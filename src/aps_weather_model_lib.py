#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
###############################################################################
# aps_weather_model_lib.py 
#
# Project:  APD TRAIN
# Purpose:  Libraries for TRAIN software
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
import glob
import saa_func_lib as saa
import numpy as np


def get_param(target):
    f = open('parms_aps.txt','r')
    for line in f:
        if target in line:
            t = line.split(':')
            return(t[1].strip())
    print "ERROR: Unable to find parameter {PARM}".format(PARM=target)
    exit(1)

def get_range(target,geo_ref_file=None):
    if geo_ref_file is not None and ".tif" in geo_ref_file:
        print "Reading {TARG} from file {FILE}".format(TARG=target,FILE=geo_ref_file)
        if target == 'region_lat_range':
            reffile = glob.glob("*wgs84.tif")[0]
            (x,y,trans,proj,data) = saa.read_gdal_file(saa.open_gdal_file(geo_ref_file))
            ullat = trans[3]
            lrlat = trans[3] + y*trans[5]
            mylist = [ullat,lrlat]
            return mylist 
        elif target == 'region_lon_range':
            reffile = glob.glob("*wgs84.tif")[0]
            (x,y,trans,proj,data) = saa.read_gdal_file(saa.open_gdal_file(geo_ref_file))
            ullon = trans[0]
            lrlon = trans[0] + x*trans[1]
            mylist = [ullon,lrlon]
            return mylist
        else:
            print "ERROR: Unknown parameter {PARM}".format(PARM=target)
            exit(1)
    else:
        print "Reading {TARG} from file params_aps.txt".format(TARG=target)
        f = open('parms_aps.txt','r')
        for line in f:
            if target in line:
                t = line.split(':')
                s = t[1].split()
                mylist = [float(s[0]),float(s[1])]
                return mylist
        print "ERROR: Unable to find parameter {PARM}".format(PARM=target)
        exit(1)

def get_date_list():
    # get list of interferograms
    datelist=[]
    for myfile in glob.glob('*_*_unw_phase.tif'):
        t = myfile.split('_')
        datelist.append(t[0])
        datelist.append(t[1])
    datelist = list(set(datelist))
    datelist.sort()
    return datelist

def times(utc,datelist):
    t_before = int(np.floor(float(utc)/21600))
    t_after = int(np.ceil(float(utc)/21600))

    temp = utc-21600*t_before
    f_after = temp/21600
    if np.isnan(f_after):
        f_after = 1
    f_before = 1 - f_after
    timelist = ['0000','0600','1200','1800','0000']
    n_SAR = len(datelist)
    date_list = []
    time_list = []
    frac_list = []
    for i in range(len(datelist)):
        date_list.append(datelist[i])
        time_list.append(timelist[t_before])
        frac_list.append(f_before)
    for i in range(len(datelist)):
        if t_after == 4:
            # must increment the date by 1 day
            t=time.strptime(date_list[i],'%Y%m%d')
            newdate=date(t.tm_year,t.tm_mon,t.tm_mday)+timedelta(1)
            thisdate = newdate.strftime('%Y%m%d')
        else:
            thisdate = date_list[i]
        date_list.append(thisdate)
        time_list.append(timelist[t_after])
        frac_list.append(f_after)
    return(date_list,time_list,frac_list)

def file_names(date,time,merra_datapath):
    output_list = []
    for i in range(len(date)):
        ymd=date[i]
        thistime = time[i]
        hour=thistime[0:2]
        path = merra_datapath + "/" + ymd + "/"
        output = "{PATH}MERRA2_{YMD}_{HOUR}.nc4".format(PATH=path,YMD=ymd,HOUR=hour)
        output_list.append(output)
    return output_list
    

