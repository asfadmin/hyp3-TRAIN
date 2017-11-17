#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
###############################################################################
# aps_merra_files.py 
#
# Project:  APD TRAIN
# Purpose:  Download MERRA2 weather data 
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
import os
from execute import execute
import argparse
import glob
from os.path import expanduser
import numpy as np
from datetime import date, timedelta
import aps_weather_model_lib as aps

def aps_merra_files(order_flag,geo_ref_file=None):  
  
    # get list of interferograms
    datelist, tmp = aps.get_date_list()
    n_files = len(datelist)

    # read URS credentials
#    home = expanduser("~")
#    myfile = home + '/' + '.merrapass'
#    f = open(myfile,'r')
#    lines = f.read()
#    t = lines.split('\n')
#    username = t[0]
#    password = t[1]
 
    # setup parameters
    utc = float(aps.get_param('UTC_sat'))
    merra_datapath = aps.get_param('merra2_datapath')

    # determine lat, lon range
    region_lat_range = aps.get_range('region_lat_range', geo_ref_file = geo_ref_file)
    region_lon_range = aps.get_range('region_lon_range', geo_ref_file = geo_ref_file)

   # add 2 degrees to lat,lon range
    S = round(min(region_lat_range)) - 2
    N = round(max(region_lat_range)) + 2
    W = round(min(region_lon_range)) - 2
    E = round(max(region_lon_range)) + 2

    # figure out which dates.times we need to download
    (date,time,frac) = aps.times(utc,datelist)

    # create download file names
    download_list = []
    for i in range(len(date)):
        ymd=date[i]
        thistime = time[i]
        hour=thistime[0:2]
        year=ymd[0:4]
        month=ymd[4:6]
        if year < 1993:
            series_type = 100
        elif year < 2001:
            series_type = 200
        elif year < 2011:
            series_type = 300
        else:
            series_type = 400
        download = "http://goldsmr5.gesdisc.eosdis.nasa.gov/daac-bin/OTF/HTTP_services.cgi?FILENAME=%2Fdata%2Fs4pa%2FMERRA2%2FM2I6NPANA.5.12.4%2F{YEAR}%2F{MONTH}%2FMERRA2_{TYPE}.inst6_3d_ana_Np.{YMD}.nc4&FORMAT=bmM0Yy8&BBOX={S}%2C{W}%2C{N}%2C{E}&TIME=1979-01-01T{HOUR}%3A00%3A00%2F1979-01-01T{HOUR}%3A00%3A00&LABEL=svc_MERRA2_{TYPE}.inst6_3d_ana_Np.{YMD}.nc4&FLAGS=&SHORTNAME=M2I6NPANA&SERVICE=SUBSET_MERRA2&LAYERS=&VERSION=1.02&VARIABLES= ".format(S=int(S),N=int(N),W=int(W),E=int(E),YEAR=year,MONTH=month,TYPE=series_type,YMD=ymd,HOUR=hour)
        download_list.append(download)
        
    # create output file names
    output_list = aps.file_names("merra2",date,time,merra_datapath)
    
    # download the data or write a file
    if int(order_flag) == 1:
        for i in range(len(date)):
            (path,myfile) = os.path.split(output_list[i])
            if not os.path.isdir(path):
                os.makedirs(path)
            if os.path.isfile(output_list[i]):
                print "File {MYFILE} exists. Skipping...".format(MYFILE=myfile)
            else:
                cmd = "wget -O{OUTFILE} \"{DOWNFILE}\"".format(OUTFILE=output_list[i],DOWNFILE=download_list[i])
                execute(cmd)
    else:
        f = open("download_files.txt","w")
        for i in range(len(date)):
            print download_list[i]
            f.write("{DOWNFILE}\n".format(DOWNFILE=download_list[i]))
        f.close()

if __name__ == '__main__':

  parser = argparse.ArgumentParser(prog='aps_merra_files',
    description='Download weather model data for a stack of interferograms')
  parser.add_argument("-g","--geo_ref_file",help="Name of file for georeferencing information")
  parser.add_argument("order_flag",help="0 - generate download file; 1 - download data")
  args = parser.parse_args()

  if args.geo_ref_file is None:
      aps_merra_files(args.order_flag)
  else:
      aps_merra_files(args.order_flag,args.geo_ref_file)
      
      

