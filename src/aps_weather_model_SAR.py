#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
###############################################################################
# aps_weather_model_SAR.py 
#
# Project:  APD TRAIN
# Purpose:  Compute SAR delays from weather model data
#           It will give the computed ZENITH hydrostatic and wet delay map in cm 
#           for the selected region.
#          
# Author:  Tom Logan
#
###############################################################################
# Copyright (c) 2017, Alaska Satellite Facility
# 
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
# 
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Library General Public License for more details.
# 
# You should have received a copy of the GNU Library General Public
# License along with this library; if not, write to the
# Free Software Foundation, Inc., 59 Temple Place - Suite 330,
# Boston, MA 02111-1307, USA.
###############################################################################

#####################
#
# Import all needed modules right away
#
#####################
import sys
import os
import re
import math
from get_dem import get_dem
from osgeo import gdal
import argparse
import aps_weather_model_lib as aps
import saa_func_lib as saa
import numpy as np
import scipy as sp
import scipy.integrate
import shutil
from aps_load_merra import aps_load_merra
from aps_weather_model_nan_check import aps_weather_model_nan_check


from joblib import Parallel, delayed
# import multiprocessing
from multiprocessing import Pool


def get_DEM(demfile,geo_ref_file=None):

    region_lat_range = aps.get_range('region_lat_range', geo_ref_file = geo_ref_file)
    region_lon_range = aps.get_range('region_lon_range', geo_ref_file = geo_ref_file)

    minlat = min(region_lat_range)
    maxlat = max(region_lat_range)
    minlon = min(region_lon_range)
    maxlon = max(region_lon_range)

    if geo_ref_file is None:
        region_res = aps.get_param('region_res')
    else:
        if ".tif" in geo_ref_file:
            print "Reading geotiff DEM file"
            (x,y,trans,proj,data) = saa.read_gdal_file(saa.open_gdal_file(geo_ref_file))
            region_res = trans[1]
        else:
            print "Unknown DEM type"
            exit(1)

    if not os.path.isfile(demfile):
        print "DEM file does not exist.  Creating it..."
        get_dem(minlon,minlat,maxlon,maxlat,demfile,0,PAP=False)
        
    (x,y,trans,proj,data) = saa.read_gdal_file(saa.open_gdal_file(demfile))
    if trans[1] != region_res:
        print "Resampling DEM file to match region_res"    
        gdal.Warp("tmp.dem",demfile,xRes=region_res,yRes=region_res,resampleAlg="cubic",dstNodata=-32767)
        shutil.move("tmp.dem",demfile)
        (x,y,trans,proj,data) = saa.read_gdal_file(saa.open_gdal_file(demfile))
    nncols = x
    nnrows = y
    demres = trans[1]
    ullat = trans[3]
    lrlat = trans[3] + y*trans[5]
    ullon = trans[0]
    lrlon = trans[0] + x*trans[1]
    xmin = min(ullon,lrlon)
    xmax = max(ullon,lrlon)
    ymin = min(ullat,lrlat)
    ymax = max(ullat,lrlat)
    
    return(data,minlon,maxlon,minlat,maxlat,demres,nncols,nnrows)

#    return (data,xmin,xmax,ymin,ymax,demres,nncols,nnrows)

def write_out_file(name,xi,yi,data):
    fid = open(name,'wb')
    data_write = np.array([np.reshape(xi,-1,order='F'),np.reshape(yi,-1,order='F'),np.reshape(data,-1,order='F')])
    data_write = np.reshape(data_write,(3,-1))
    data_write = np.transpose(data_write)
    data_write.tofile(fid)
    fid.close()

def inter2d(xivec,yivec,lonlist_matrix,latlist_matrix,cdstack,n):
    tck = sp.interpolate.bisplrep(lonlist_matrix,latlist_matrix,cdstack[:,:,n],s=3)
    newslice = sp.interpolate.bisplev(yivec,xivec,tck)
    return newslice
#    cdstack_interp_dry[:,:,n]=np.flipud(newslice)

def aps_weather_model_SAR(demfile=None,geo_ref_file=None):

    # Constants
    # Parameters from Jolivet et al 2011, GRL

    Rd = 287.05            # [j/kg/K] specific gas constants dry air
    Rv = 461.495           # [j/kg/K] water vapour
    k1 = 0.776             # [K/Pa]
    k2 = 0.716             # [K/Pa]
    k3 = 3.75e3            # [K^2/Pa]
    g0 = 9.80665
    Rmax = 6378137 
    Rmin = 6356752

#    print "k1 {}".format(k1)
#    print "k2 {}".format(k2)
#    print "k3 {}".format(k3)
    
    zref = 15000       # zref for integration- tropopause
    zincr = 15         # z increment spacing for vertical interpolation before integration
    vertres = 100      # vertical resolution for comparing dem to vert profiles of delay

    XI= [int(i) for i in range(0,zref+1,zincr)]
#    print "len(XI) {}".format(len(XI))
    
    path = aps.get_param('merra_datapath')
    demfile = aps.get_param('demfile')
    utc = float(aps.get_param('UTC_sat'))
    dates, tmp = aps.get_date_list()
    n_dates = len(dates)
     
    dem,xmin,xmax,ymin,ymax,smpres,nncols,nnrows = get_DEM(demfile,geo_ref_file)
    
#    print xmin, xmax
#    print ymin, ymax
#    print smpres
#    print nncols, nnrows
#    print dem[:30]
    
    lonmin = np.floor(xmin)-1
    lonmax= np.ceil(xmax)+1
    latmin = np.floor(ymin)-1 
    latmax = np.ceil(ymax)+1

    maxdem = np.ceil(np.max(np.max(dem))/100)*100+100
    cdslices = int(maxdem/vertres) + 1
    cdI = [x for x in range(0,int(maxdem)+1,vertres)]
#    print "cdI {}".format(cdI)
    
    rounddem = np.round(dem/vertres)
    rounddem[dem<0] = 0
    rounddem[np.isnan(dem)] = 0

    if cdslices != len(cdI):
       print "ERROR: length mistmatch"
       print "cdslices {}".format(cdslices)
       print "len(cdI) {}".format(len(cdI))
       print "cdI {}".format(cdI)
       exit(1)

    print "Interpolate to a maximum dem height of {HEIGHT}".format(HEIGHT=maxdem)
    print "Using {} slices".format(cdslices)
    
    date,time,frac = aps.times(utc,dates)
   
    input_file_names = aps.file_names(date,time,path)
    length = len(input_file_names)/2
    
    for d in range(length):
        no_data = 0
        for kk in range(2):
            if kk==0:
                wfile = input_file_names[d]
            else:
                wfile = input_file_names[d+length]
                
            if not os.path.isfile(wfile):
                no_data = no_data + 1
            else:
                print "Loading input weather file: {}".format(wfile)
            
                # Load the weather model data
                (Temp,e,H,P,longrid,latgrid,xx,yy,lon0360_flag) = aps_load_merra(wfile) 
                
#                print "H[:,:,0] {}".format(H[:,:,0])
#                print "H[:,:,1] {}".format(H[:,:,1])
#                print "H[:,:,2] {}".format(H[:,:,2])
#                print "H[:,:,3] {}".format(H[:,:,3])
#                print "H[:,:,4] {}".format(H[:,:,4])
#                print "H[:,:,5] {}".format(H[:,:,5])
    

#                print "H {}".format(H)
#                print "P {}".format(P)
                
 		# deal with NANs
                (Temp,e) = aps_weather_model_nan_check(Temp,e,P,longrid,latgrid)                

#                print "Temp[:,:,0] {}".format(Temp[:,:,0])
#                print "Temp[:,:,1] {}".format(Temp[:,:,1])
#                print "Temp[:,:,2] {}".format(Temp[:,:,2])
#                print "Temp[:,:,3] {}".format(Temp[:,:,3])
#                print "Temp[:,:,4] {}".format(Temp[:,:,4])
#                print "Temp[:,:,5] {}".format(Temp[:,:,5])

#                print "e[:,:,0] {}".format(e[:,:,0])
#                print "e[:,:,1] {}".format(e[:,:,1])
#                print "e[:,:,2] {}".format(e[:,:,2])
#                print "e[:,:,3] {}".format(e[:,:,3])
#                print "e[:,:,4] {}".format(e[:,:,4])
#                print "e[:,:,5] {}".format(e[:,:,5])


		# define weather model grid nodes
                latlist = np.reshape(latgrid[:,:,1],-1)
                lonlist = np.reshape(longrid[:,:,1],-1)
                
#                print latlist
#                print lonlist
                
#                print "latlist shape {}".format(latlist.shape)
#                print "lonlist shape {}".format(lonlist.shape)
                
                xlist = np.reshape(xx,-1,order='F')
                ylist = np.reshape(yy,-1,order='F')
                
#                print xlist
#                print ylist
                
                lat_res = np.abs(np.diff(np.unique(latgrid)))*1.5
                lat_res = lat_res[0]
                lon_res = np.abs(np.diff(np.unique(longrid)))*1.5
                lon_res = lon_res[0]

#                print "lat res {}".format(lat_res)
#                print "lon res {}".format(lon_res)
                
                if lon0360_flag == 'y':
                    if xmin < 0:
                        xmin = xmin + 360
                    if xmax < 0:
                        xmax = xmax + 360
                
                # ix = np.where(ymin-lat_res<= latlist and latlist <= ymax+lat_res and xmin-lon_res<= lonlist and lonlist<= xmax+lon_res)
                ix1 = np.where(ymin-lat_res<= latlist)
                ix2 = np.where(latlist <= ymax+lat_res)
                ix3 = np.intersect1d(ix1,ix2)
                ix1 = np.where(xmin-lon_res<= lonlist)
                ix2 = np.where(lonlist<= xmax+lon_res)
                ix4 = np.intersect1d(ix1,ix2)
                ix = np.intersect1d(ix3,ix4)
 
#                print "lat_res {}".format(lat_res)
#                print "lon_res {}".format(lon_res)
#                print "ymin    {}".format(ymin)
#                print "ymax    {}".format(ymax)
#                print "xmin    {}".format(xmin)
#                print "xmax    {}".format(xmax)
#                print ix
                
                xlist = xlist[ix]
                ylist = ylist[ix]
                latlist = latlist[ix]
                lonlist = lonlist[ix]
                ulatlist = np.unique(latlist)
                ulonlist = np.unique(lonlist)
                numy = int(len(ulatlist))
                numx = int(len(ulonlist))
                uxlist = np.unique(xlist)
                uylist = np.unique(ylist)
  
#  		print xlist
#                print ylist
#                print latlist
#                print lonlist
#                print ulatlist
#                print ulonlist
#                print numx, numy
#                print uxlist
#                print uylist
  
                # map of g with latitude
                g = 9.80616*(1. - 0.002637*np.cos(2*np.deg2rad(latgrid)) + 0.0000059* np.square(np.cos(2.*np.deg2rad(latgrid))))
                
#                print g[:,:,0]
                
                
                # map of Re with latitude
                Re = np.sqrt(1./((np.square(np.cos(np.deg2rad(latgrid)))/np.square(Rmax)) + (np.square(np.sin(np.deg2rad(latgrid)))/np.square(Rmin))))
                
#                print Re[:,:,0]
                
                # Calculate Geometric Height, Z
                Z = (H*Re)/(g/g0*Re - H)

#                print "HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH"
#                print H[:,:,0]
#                print Z[:,:,0]

                midx = int(round(np.mean(uxlist)))
                midy = int(round(np.mean(uylist)))
                
#                print midx,midy
                
                glocal = g[midx,midy,0]
                Rlocal = Re[midx,midy,0]

#                print "glocal {}".format(glocal)
#                print "Rlocal {}".format(Rlocal)
                
                cdstack = np.zeros((numy,numx,cdslices))
                cdstack_dry = np.zeros((numy,numx,cdslices))
                cdstack_wet = np.zeros((numy,numx,cdslices))

                # Interpolate Temp P and e from 0:20:15000 m
                # then integrate using trapz to estimate delay as function of height
                for i in range(numx):
                    for j in range(numy):
                        xn = uxlist[i]
                        yn = uylist[j]
                        
#                        print "xn {}".format(xn)
#                        print "yn {}".format(yn)
                        
                        #interpolate at res zincr before doing integration
                        X = np.squeeze(Z[xn,yn,:])
#                        print X
                        
                        Ye = np.squeeze(e[xn,yn,:])
#                        print Ye
                        
                        f = sp.interpolate.splrep(X,Ye)
                        YeI = sp.interpolate.splev(XI,f,der=0)
                        YeI = YeI * 100
#                        print YeI[:30]
                        
                        Yp = np.squeeze(P[yn,xn,:])
                        f = sp.interpolate.splrep(X,Yp)
                        YPI = sp.interpolate.splev(XI,f,der=0)
                        YPI = YPI * 100
#                        print YPI[:30]


                        YT = np.squeeze(Temp[yn,xn,:])
                        f = sp.interpolate.splrep(X,YT)
                        YTI = sp.interpolate.splev(XI,f,der=0)
#                        print YTI[:30]

                        tmp1 = (k2-(Rd*k1/Rv))*YeI/YTI + k3*YeI/(YTI*YTI)
#                        print "tmp1 {}".format(tmp1[:30])
                        Lw = (10**-6)*-1*np.flipud(sp.integrate.cumtrapz(np.flipud(tmp1),np.flipud(XI),initial=0))
#                        print "Lw {}".format(Lw[:30])
                        
                        f = sp.interpolate.splrep(XI,Lw)
                        LwI = sp.interpolate.splev(cdI,f,der=0)
			
#                        print "XI {}".format(XI[:30])
#                        print "cdI {}".format(cdI[:30])
                        
                        Ld = 10**-6*((k1*Rd/glocal)*(YPI-YPI[zref/zincr]))
                        f = sp.interpolate.splrep(XI,Ld)
                        LdI = sp.interpolate.splev(cdI,f,der=0)

#			print("LwI {}".format(LwI))
                        
#                        print LdI
#                        print LwI
                        
                        cdstack_dry[j,i,:] = LdI
                        cdstack_wet[j,i,:] = LwI
                

                xsmpres = (xmax-xmin)/nncols
                ysmpres = (ymax-ymin)/nnrows
                
                xivec = np.linspace(xmin+0.5*xsmpres,xmax-0.5*xsmpres,nncols)
                yivec = np.linspace(ymin+0.5*ysmpres,ymax-0.5*ysmpres,nnrows)
	
#                print "nncols {}".format(nncols)
#                print "nnrows {}".format(nnrows)
        
#                print "Length of xivec {}".format(len(xivec))
#                print "Length of yivec {}".format(len(yivec))                

                [xi,yi] = np.meshgrid(xivec,yivec)
                ix_temp = np.diff(lonlist)
                ix_temp = np.where(ix_temp!=0)
                ixi_temp = ix_temp[0] + 1
                lonlist_matrix = np.reshape(lonlist,(ixi_temp[0],-1),order='F')
                latlist_matrix = np.reshape(latlist,(ixi_temp[0],-1),order='F')

#                print lonlist_matrix
#                print latlist_matrix
                
                cdstack_interp_dry = np.zeros((nnrows,nncols,cdslices))
                sys.stdout.write("processing dry stack")
                sys.stdout.flush()
                for n in range(cdslices):
#                    f = scipy.interpolate.interp2d(lonlist_matrix,latlist_matrix,cdstack_dry[:,:,n])
#                    newslice = f(xivec,yivec)

                    tck = sp.interpolate.bisplrep(lonlist_matrix,latlist_matrix,cdstack_dry[:,:,n],s=3)
                    newslice = sp.interpolate.bisplev(yivec,xivec,tck)

                    cdstack_interp_dry[:,:,n]=np.flipud(newslice)
                    sys.stdout.write(".")
                    sys.stdout.flush()
                print


#                Parallel(n_jobs=8)(delayed(inter2d)(xivec,yivec,lonlist_matrix,latlist_matrix,cdstack,cdstack_interp_dry,n)
#                           for n in range(cdslices))

#                pool = Pool(processes=8)
#                cdstack_interp_dry = [pool.apply(inter2d,args=(xivec,yivec,lonlist_matrix,latlist_matrix,cdstack,n)) for n in range(cdslices)]
#                print "cdstack_interp_dry shape {}".format(cdstack_interp_dry.shape())


#                print cdstack_dry[:,:,0]
#                print cdstack_interp_dry[:30,0,0]

                # del cd_stack_interp_wet
                cdstack_interp_wet = np.zeros((nnrows,nncols,cdslices))
                sys.stdout.write("processing wet stack")
                sys.stdout.flush()
                
                lonlist_matrix = lonlist_matrix[0,:]
                latlist_matrix = latlist_matrix[:,0]
                
#                print lonlist_matrix
#                print latlist_matrix
                
                for n in range(cdslices):
#                    f = sp.interpolate.interp2d(lonlist_matrix,latlist_matrix,cdstack_wet[:,:,n])

                    f = sp.interpolate.RectBivariateSpline(latlist_matrix,lonlist_matrix,cdstack_wet[:,:,n],kx=1,ky=1,s=0)
                    newslice = f(yivec,xivec)
                    
#                    print newslice

#                    tck = sp.interpolate.bisplrep(lonlist_matrix,latlist_matrix,cdstack_wet[:,:,n],s=5)
#                    newslice = sp.interpolate.bisplev(yivec,xivec,tck)

                    cdstack_interp_wet[:,:,n]=np.flipud(newslice)
                    sys.stdout.write(".")
                    sys.stdout.flush()
                print

#                print cdstack_wet[:,:,0]
#                print cdstack_interp_wet[:30,0,0]
                
#                cdstack_wet[:,:,0].tofile("cstack_wet.bin")
#                cdstack_interp_wet[:,:,0].tofile("cdstack_interp_wet.bin")
                
                
                # keeping the coordinates in the same grid as the data
                xi = np.flipud(xi)
                yi = np.flipud(yi)
                
                # Pull out delays from cdstack layers that match dem heights
                wetcorrection = np.ones((nnrows,nncols))
                hydrcorrection = np.ones((nnrows,nncols))
                
#                print dem[:30,0]
#                print rounddem[:30,0]
                
                for i in range(nnrows):
                    for j in range(nncols):
                        wetcorrection[i,j] = cdstack_interp_wet[i,j,int(rounddem[i,j])]
                        hydrcorrection[i,j] = cdstack_interp_dry[i,j,int(rounddem[i,j])]
                
#                print wetcorrection[:30,0]
                
#                del cdstack_interp_wet
#                del cdstack_interp_dry
                
                if kk == 0:
                    wetcorrection1 = wetcorrection
                    drycorrection1 = hydrcorrection
                elif kk == 1:
                    wetcorrection2 = wetcorrection
                    drycorrection2 = hydrcorrection
                del wetcorrection 
                del hydrcorrection
        
        minf = wetcorrection2.min()
        maxf = wetcorrection2.max()
        
        where = np.argmin(wetcorrection2)
        
#        print "minimum wetcorrection2 {MINF} at {X}".format(MINF=minf,X=where)
#        print "maximum wetcorrection2 {}".format(maxf)
        
        
        if no_data == 0:
            
            print "Writing output files"
            
            outfile_wet_before = path + "/" + date[d] + "/" + date[d] + '_ZWD_before.xyz'
            outfile_wet_after = path + "/" + date[d] + "/" + date[d] + '_ZWD_after.xyz'
            
            outfile_dry_before = path + "/" + date[d] + "/" + date[d] + '_ZHD_before.xyz'
            outfile_dry_after = path + "/" + date[d] + "/" + date[d] + '_ZHD_after.xyz'            
            
            write_out_file(outfile_wet_before,xi,yi,wetcorrection1)
            write_out_file(outfile_wet_after,xi,yi,wetcorrection2)
            
            write_out_file(outfile_dry_before,xi,yi,drycorrection1)
            write_out_file(outfile_dry_after,xi,yi,drycorrection2)
                       
            # Output wet correction
            wetcorrection = wetcorrection1*frac[d]+wetcorrection2*frac[d+length]
            del wetcorrection1 
            del wetcorrection2
            wetcorrection = wetcorrection * 100     # delay in CM
            outfile = path + "/" + date[d] + "/" + date[d] + '_ZWD.xyz'
            write_out_file(outfile,xi,yi,wetcorrection)
            del wetcorrection
            
            # Output hydrostatic correction
            hydrcorrection = drycorrection1*frac[d]+drycorrection2*frac[d+length]
            del drycorrection1 
            del drycorrection2
            hydrcorrection = hydrcorrection * 100   # delay in CM
            hydroutfile = path + "/" + date[d] + "/" + date[d] + '_ZHD.xyz'
            write_out_file(hydroutfile,xi,yi,hydrcorrection)
            del hydrcorrection
            
            print "{D} completed out of {L}".format(D=d+1,L=length)
        else:
            print "{D} completed out of {L} (NO DATA)".format(D=d+1,L=length)
        
    
if __name__ == '__main__':

  parser = argparse.ArgumentParser(prog='aps_wetaher_model_SAR',
    description='Calcuate zenith wet and hydrostatic delays')
  parser.add_argument("-d","--dem",help="Name of DEM file to use (geotiff)")
  parser.add_argument("-g","--geo_ref_file",help="Name of file for georeferencing information")

  args = parser.parse_args()
  aps_weather_model_SAR(demfile=args.dem,geo_ref_file=args.geo_ref_file)

