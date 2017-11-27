#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
###############################################################################
# aps_weather_model_lib.py 
#
# Project:  APD TRAIN
# Purpose:  Libraries for TRAIN software
#          
# Author:  Tom Logan, adapted from David Bekaert's matlab code
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
import os
import re
import saa_func_lib as saa
import numpy as np
import datetime
import time

def timestamp(date):
    return time.mktime(date.timetuple())

def get_param(target):
    f = open('parms_aps.txt','r')
    print "Reading {TARG} from file parms_aps.txt".format(TARG=target)
    for line in f:
        if target in line:
            t = line.split(':')
            return(t[1].strip())
    print "ERROR: aps_weather_model was unable to find parameter {PARM} in the parms_aps.txt file".format(PARM=target)
    print "ERROR: Please add this parameter to your parms_aps.txt file and try again"
    exit(1)

def get_range(target,geo_ref_file=None):
    if geo_ref_file is not None and ".tif" in geo_ref_file:
        print "Reading {TARG} from file {FILE}".format(TARG=target,FILE=geo_ref_file)
        if target == 'region_lat_range':
            (x,y,trans,proj,data) = saa.read_gdal_file(saa.open_gdal_file(geo_ref_file))
            ullat = trans[3]
            lrlat = trans[3] + y*trans[5]
            mylist = [ullat,lrlat]
            return mylist 
        elif target == 'region_lon_range':
            (x,y,trans,proj,data) = saa.read_gdal_file(saa.open_gdal_file(geo_ref_file))
            ullon = trans[0]
            lrlon = trans[0] + x*trans[1]
            mylist = [ullon,lrlon]
            return mylist
        else:
            print "ERROR: Unknown parameter {PARM} in get_range function".format(PARM=target)
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
    origin = get_param('date_origin')
    if origin == 'asf':
        for myfile in glob.glob('*_*_unw_phase.tif'):
            t = myfile.split('_')
            datelist.append(t[0])
            datelist.append(t[1])
    elif origin == 'file':
        fname = get_param('ifgday_file')
        if os.path.isfile(fname):
            f = open(fname,'r')
            for line in f:
                t = re.split(' ',line)
                datelist.append(t[0].strip())
                datelist.append(t[1].strip())
            f.close()
        else:
            print "ERROR: Unable to find ifgday file {}".format(fname)
            exit(1)
    else:
        print "ERROR: Unknown date_origin type ({}) read from parms_aps.txt file".format(origin)
        exit(1)

    shortlist = list(set(datelist))
    shortlist.sort()
    if len(shortlist) == 0:
        print "ERROR: No dates found; Nothing to do"
        print ""
        if origin == 'asf':
            print "No files of the form *_unw_phase.tif were found in {}".format(os.getcwd())
            print "You must place your unwrapped phase geotiffs in the directory where you ran the code."
        exit(1)

    return shortlist, datelist

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

def file_names(model_type,date,time,datapath):
    output_list = []
    for i in range(len(date)):
        ymd=date[i]
        thistime = time[i]
        hour=thistime[0:2]
        path = datapath + "/" + ymd + "/"
        output = "{PATH}{TYPE}_{YMD}_{HOUR}.nc4".format(PATH=path,TYPE=model_type.upper(),YMD=ymd,HOUR=hour)
        output_list.append(output)
    return output_list
    

