#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
###############################################################################
# aps_weather_model.py 
#
# Project:  APD TRAIN
# Purpose:  Compute SAR delays from a stack of interferograms.
#           This is the main driver program for the APS weather correction.
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
import argparse
import aps_weather_model_lib as aps
from aps_merra_files import aps_merra_files
from aps_weather_model_SAR import aps_weather_model_SAR
from aps_weather_model_INSAR import aps_weather_model_INSAR
from aps_weather_model_diff import aps_weather_model_diff
import datetime
import time

def aps_weather_model(model_type,start,end,geo_ref_file=None):
    start_time = datetime.datetime.now()

    if start == 0:
        print "Step 0: Determine files to download"
        if model_type=='era':
            print "ERROR: ERA is still under construction"
            exit(1)
        elif model_type=='merra':
            print "ERROR: MERRA is still under construction"
            exit(1)
        elif model_type == 'merra2':
            aps_merra_files(0,geo_ref_file)
        else:
            print "ERROR: Unknown model type {}".format(model_type)
            exit(1)
    if start <= 1 and end >= 1:
        print "Step 1: Order and download {type} data files".format(type=model_type)
        if model_type=='era':
            print "ERROR: ERA is still under construction"
            exit(1)
        elif model_type=='merra':
            print "ERROR: MERRA is still under construction"
            exit(1)
        elif model_type == 'merra2':
            aps_merra_files(1,geo_ref_file)
        else:
            print "ERROR: Unknown model type {}".format(model_type)
            exit(1)
    if start <= 2 and end >= 2:
        print "Step 2: Compute the {type} tropospheric (zenith) delay for individual dates".format(type=model_type)
        aps_weather_model_SAR(model_type,geo_ref_file)
    if start <= 3 and end >= 3:
        print "Step 3: Computes {type} tropospheric (slant) delay for inteferograms".format(type=model_type)
        aps_weather_model_INSAR(model_type,geo_ref_file)
    if start <= 4 and end >= 4:
        print "Step 4: Subtract delay from interferograms"
        aps_weather_model_diff()
    
    print "Done!!!"
    print "Total processing time {}".format(aps.timestamp(datetime.datetime.now())-aps.timestamp(start_time))
    
if __name__ == '__main__':

  parser = argparse.ArgumentParser(prog='aps_weather_model',
    description="Calcuate tropospheric (slant) delay for a stack of interferograms: \n\
                 0 - Set up list of weather data files to download; \n\
                 1 - Download weather data; \n\
                 2 - Calculate wet and hydrostatic zenith delays; \n\
                 3 - Calculate the SAR delays; \n\
                 4 - Subtract calculated delay from interferograms")
  parser.add_argument("sstep",help="Start step to begin processing (0-4)")
  parser.add_argument("estep",help="End step to stop processing (0-4)")
  parser.add_argument("-g","--geo_ref_file",help="Name of file for georeferencing information")
  args = parser.parse_args()
  
  aps_weather_model("merra2",int(args.sstep),int(args.estep),geo_ref_file=args.geo_ref_file)
