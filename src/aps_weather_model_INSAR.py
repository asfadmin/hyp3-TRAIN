#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
###############################################################################
# aps_weather_model_INSAR.py 
#
# Project:  APD TRAIN
# Purpose:  Convert zenith hydrostatic and wet delays into phase corrections
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
import os
import argparse
import aps_weather_model_lib as aps
import saa_func_lib as saa
import numpy as np
import scipy as sp
import scipy.integrate
import scipy.interpolate
import datetime
import time
from osgeo import gdal

def get_ll_mat(geo_ref_file):
    if geo_ref_file is not None:
        (x,y,trans,proj,data) = saa.read_gdal_file(saa.open_gdal_file(geo_ref_file))
        lonlat = np.zeros((x*y,2))
        ullat = trans[3]
        ullon = trans[0]
	lon = np.zeros(x)
        for i in range(x):
            lon[i] = ullon + i*trans[1]
        lat = np.zeros(y)
        for i in range(y):
            lat[i] = ullat + i*trans[5]
	[xi,yi] = np.meshgrid(lon,lat)
        xi = np.reshape(xi,-1)
        yi = np.reshape(yi,-1)
	lonlat[:,0] = xi
        lonlat[:,1] = yi
    else:
        print "llmat under construction; must use geo_ref_file for now"
        exit(1)
    return lonlat

def load_weather_model_SAR(infile, output_grid, chunk_flag=True, chunk_size=0.5, geo_ref_file=None):
    print "Reading file {}".format(infile)
    fid = open(infile,'rb')
    data = np.fromfile(fid)
    data = np.reshape(data, (-1,3))
       
    if chunk_flag:
        xy_min0 = np.floor(np.min(output_grid[:,0]))
        xy_min1 = np.floor(np.min(output_grid[:,1]))
        xy_max0 = np.ceil(np.max(output_grid[:,0]))
        xy_max1 = np.ceil(np.max(output_grid[:,1]))
        lon_grid = np.arange(xy_min0,xy_max0-chunk_size+0.1,chunk_size)
        lat_grid = np.arange(xy_min1,xy_max1-chunk_size+0.1,chunk_size)
        lon_grid,lat_grid = np.meshgrid(lon_grid,lat_grid)
        lon_grid = np.ndarray.flatten(lon_grid,order='F')
        lat_grid = np.ndarray.flatten(lat_grid,order='F')
        z_output = np.zeros(len(output_grid[:,0]))      
        print "    Buffering in {} deg chunks...".format(chunk_size)
        for k in range(len(lat_grid)):
            ix_temp1 = np.where(output_grid[:,0]>=lon_grid[k])
            ix_temp2 = np.where(output_grid[:,0]<=lon_grid[k]+chunk_size)
            ix_temp3 = np.where(output_grid[:,1]>=lat_grid[k])
            ix_temp4 = np.where(output_grid[:,1]<=lat_grid[k]+chunk_size)
            ix_temp5 = np.intersect1d(ix_temp1,ix_temp2)
            ix_temp6 = np.intersect1d(ix_temp3,ix_temp4)
            ix_temp = np.intersect1d(ix_temp5,ix_temp6)
            
            if len(ix_temp) > 0:
                ix_temp1 = np.where(data[:,0]>=lon_grid[k]-chunk_size/2)
                ix_temp2 = np.where(data[:,0]<=lon_grid[k]+chunk_size+chunk_size/2)
                ix_temp3 = np.where(data[:,1]>=lat_grid[k]-chunk_size/2)
                ix_temp4 = np.where(data[:,1]<=lat_grid[k]+chunk_size+chunk_size/2)
                ix_temp5 = np.intersect1d(ix_temp1,ix_temp2)
                ix_temp6 = np.intersect1d(ix_temp3,ix_temp4)
                ix_temp_full = np.intersect1d(ix_temp5,ix_temp6)
                
                input_grid = np.zeros((len(ix_temp_full),2))
                input_grid[:,0] = data[ix_temp_full,0]
                input_grid[:,1] = data[ix_temp_full,1]
                
                tmp_output = sp.interpolate.griddata(input_grid,data[ix_temp_full,2],output_grid[ix_temp])
                
                if len(ix_temp) == len(tmp_output):
                    # print "storing data in z_output"
                    z_output[ix_temp] = tmp_output
            print "        {K}/{TOT}".format(K=k+1,TOT=len(lat_grid))
    else:
        input_grid = np.zeros((len(data[:,0]),2))
        input_grid[:,0] = data[:,0]
        input_grid[:,1] = data[:,1]
        z_output = sp.interpolate.griddata(input_grid,data[:,2],output_grid)
    return z_output


def load_weather_model_SAR2(infile,geo_ref_file):
    demfile = aps.get_param('demfile')
    (x,y,trans,proj,data) = saa.read_gdal_file(saa.open_gdal_file(geo_ref_file))
    (x1,y1,trans1,proj1,data1) = saa.read_gdal_file(saa.open_gdal_file(demfile))
    print "Reading file {}".format(infile)
    fid = open(infile,'rb')
    data = np.fromfile(fid)
    data = np.reshape(data, (-1,3))
    
    if len(data[:,2]) != x1*y1:
        print "Error: data size mismatch"
        exit(1)
    data = data[:,2]
    data = np.reshape(data,(x1,y1))
    data = np.transpose(data)
    outfile = os.path.basename(infile.replace(".xyz",".tif"))
    saa.write_gdal_file_float(outfile,trans1,proj1,data)
    outfile2 = "regrid_" + outfile
    
    ullat = trans[3]
    lrlat = trans[3] + y*trans[5]
    ullon = trans[0]
    lrlon = trans[0] + x*trans[1]
    output_bounds = [min(lrlon,ullon),min(lrlat,ullat),max(lrlon,ullon),max(lrlat,ullat)]
    gdal.Warp(outfile2,outfile,width=x,height=y,outputBounds=output_bounds,resampleAlg='bilinear')
    (x,y,trans,proj,data) = saa.read_gdal_file(saa.open_gdal_file(outfile2))
    
    os.remove(outfile)
    os.remove(outfile2)
    
    data = data.flatten()
    return data

def aps_weather_model_INSAR(geo_ref_file=None):

    start = datetime.datetime.now()
    print "Calculating phase delays"
    
    path = aps.get_param('merra_datapath')
    dates, datelist = aps.get_date_list()
    n_dates = len(dates)
    n_igrams = len(datelist)/2
    dates_master = datelist[0:len(datelist)+1:2]
    dates_slave = datelist[1:len(datelist)+1:2]
    idx_master = []
    idx_slave = []
    for i in range(n_igrams):
         idx_master.append(dates.index(dates_master[i]))
         idx_slave.append(dates.index(dates_slave[i]))

    lamb = aps.get_param('lambda')
    lamb = float(lamb)*100
    look_angle = float(aps.get_param('look_angle'))
    lonlat = get_ll_mat(geo_ref_file)
    
    d_wet = np.zeros((len(lonlat[:,0]),n_dates))
    d_wet[:] = np.nan
    d_hydro = np.zeros((len(lonlat[:,0]),n_dates))
    d_hydro[:] = np.nan    
    counter = 0
    ix_no_weather_model_data = []
    
    table = np.zeros(n_dates)
    for k in range(n_dates):
        if os.path.isfile("hydro_correction{}.bin".format(k)) and os.path.isfile("wet_correction{}.bin".format(k)):
            table[k] = 1
    if np.sum(table)!=n_dates:
        for k in range(n_dates):
            model_file_wet = path + "/" + dates[k] + "/" + dates[k] + '_ZWD.xyz'
            model_file_hydro = path + "/" + dates[k] + "/" + dates[k] + '_ZHD.xyz'
            
            if os.path.isfile(model_file_wet) and os.path.isfile(model_file_hydro):
                lasttime = datetime.datetime.now()
                # d_hydro[:,k] = load_weather_model_SAR(model_file_hydro,lonlat,geo_ref_file=geo_ref_file)
                
                d_hydro[:,k] = load_weather_model_SAR2(model_file_hydro,geo_ref_file)
                
#                outfile = "hydro_correction" + str(k) + ".bin"
#                d_hydro[:,k].tofile(outfile)
                print("Time to read and interpolate file {}".format(aps.timestamp(datetime.datetime.now())-aps.timestamp(lasttime)))
                lasttime = datetime.datetime.now()
                # d_wet[:,k] = load_weather_model_SAR(model_file_wet,lonlat,geo_ref_file=geo_ref_file)
                
                d_wet[:,k] = load_weather_model_SAR2(model_file_wet,geo_ref_file)
                
#                outfile = "wet_correction" + str(k) + ".bin"
#                d_wet[:,k].tofile(outfile)
                print("Time to read and interpolate file {}".format(aps.timestamp(datetime.datetime.now())-aps.timestamp(lasttime)))
                counter = counter + 1
            else:
                print "no correction for {}. Exiting.".format(date[k])
                ix_no_weather_model_data.append(k)
              
        print "{CNT} out of {NDATES} SAR images have a tropospheric delay estimated".format(CNT=counter,NDATES=n_dates)

    else:
        print "Reading in previously calculated intermediate products"
        for k in range(n_dates):
            lasttime = datetime.datetime.now()
            infile = "hydro_correction" + str(k) + ".bin"
            d_hydro[:,k] = np.fromfile(infile,dtype='double')
            print("Time to read file {}".format(aps.timestamp(datetime.datetime.now())-aps.timestamp(lasttime)))
            lasttime = datetime.datetime.now()
            infile = "wet_correction" + str(k) + ".bin"
            d_wet[:,k] = np.fromfile(infile,dtype='double')
            print("Time to read file {}".format(aps.timestamp(datetime.datetime.now())-aps.timestamp(lasttime)))
            
    d_total = d_hydro + d_wet
    
    d_total = d_total / np.cos(look_angle)
    d_hydro = d_hydro / np.cos(look_angle)
    d_wet = d_wet / np.cos(look_angle)
    
    ph_SAR = -4 * np.pi/lamb*d_total
    ph_SAR_hydro = -4 * np.pi/lamb*d_hydro
    ph_SAR_wet = -4 * np.pi/lamb*d_wet
    
    del d_total
    del d_hydro
    del d_wet
    
    ph_tropo = np.zeros(len(lonlat[:,0]))
    ph_tropo_hydro = np.zeros(len(lonlat[:,0]))
    ph_tropo_wet = np.zeros(len(lonlat[:,0]))
    
    for k in range(n_igrams):
        ph_tropo[:] = ph_SAR[:,idx_master[k]] - ph_SAR[:,idx_slave[k]]
        ph_tropo_hydro[:] = ph_SAR_hydro[:,idx_master[k]] - ph_SAR_hydro[:,idx_slave[k]]
        ph_tropo_wet[:] = ph_SAR_wet[:,idx_master[k]] - ph_SAR_wet[:,idx_slave[k]]
        
        outfile = dates_master[k] + "_" + dates_slave[k] + "_correction.bin"
        ph_tropo.tofile(outfile)

        outfile = dates_master[k] + "_" + dates_slave[k] + "_hydro_correction.bin"
        ph_tropo_hydro.tofile(outfile)
        
        outfile = dates_master[k] + "_" + dates_slave[k] + "_wet_correction.bin"
        ph_tropo_wet.tofile(outfile)

    print "Total processing time {}".format(aps.timestamp(datetime.datetime.now())-aps.timestamp(start))

if __name__ == '__main__':

  parser = argparse.ArgumentParser(prog='aps_weather_model_INSAR',
    description='Convert zenith delay to phase delay')
  parser.add_argument("-g","--geo_ref_file",help="Name of file for georeferencing information")

  args = parser.parse_args()
  aps_weather_model_INSAR(geo_ref_file=args.geo_ref_file)

