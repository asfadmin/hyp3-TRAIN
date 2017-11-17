#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
###############################################################################
# aps_weather_model_diff.py 
#
# Project:  APD TRAIN
# Purpose:  Subtract phase delays from interferograms
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
import re
import argparse
from osgeo import gdal
import glob
import numpy as np
import saa_func_lib as saa

def aps_weather_model_diff():
    
    for myfile in glob.glob('*_*_unw_phase.tif'):
        
        # read interferogram phase file
        print "Reading file {}".format(myfile)
        (x,y,trans,proj,data) = saa.read_gdal_file(saa.open_gdal_file(myfile))
        flat_data = data.flatten()
        
        # get dates of this interferogram
        t = re.split("_",myfile)
        date1 = t[0]
        date2 = t[1]
        
	print "found dates {a} {b}".format(a=date1,b=date2)

        # read correction file
        name = date1 + "_" + date2 + "_correction.bin"
        print "Reading file {}".format(name)
        correction = np.fromfile(name)
        np.putmask(correction,data==0,0)
                
        # subtract
        out = flat_data - correction
        outdata = np.reshape(out,(y,x))
        
        # write dfiference file
        name = date1 + "_" + date2 + "_unw_phase_corrected.tif"
        print "Writing file {}".format(name)
        saa.write_gdal_file_float(name,trans,proj,outdata)



if __name__ == '__main__':

  parser = argparse.ArgumentParser(prog='aps_weather_model_diff',
    description='Subtract phase delay from interferogram phase')

  args = parser.parse_args()
  aps_weather_model_diff()


