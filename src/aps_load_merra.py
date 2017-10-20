#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
###############################################################################
# aps_load_merra.py 
#
# Project:  APD TRAIN
# Purpose:  Load MERRA2 netcdf file
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
import re
import math
from execute import execute
from osgeo import gdal
import argparse
import glob
import aps_weather_model_lib as aps
import saa_func_lib as saa
import shutil
from netCDF4 import Dataset as NetCDFFile 
import matplotlib.pyplot as plt
import numpy as np

# from mpl_toolkits.basemap import Basemap


def aps_load_merra(wfile):

    # Set a couple of constants
    missing_value = 999999986991104
    Rv= 461.495           # [j/kg/K] water vapour
    Rd= 287.05            # [j/kg/K] specific gas constants dry air

    # Read netCDF file
    nc = NetCDFFile(wfile)
    qv = nc.variables['QV'][:]
    H = nc.variables['H'][:]
#    g0 = 9.80665
#    Geopot = H*g0
    lons = nc.variables['lon'][:]
    lats = nc.variables['lat'][:]
    Temp = nc.variables['T'][:]
    Plevs = nc.variables['lev'][:]
    time = nc.variables['time'][:]
    
    # Convert time to hours
    time = time / 60
    
#    print time
#    print lons
#    print lats
#    print Temp
    
    # Update no data values with NANs
    Temp[Temp==missing_value]=np.nan
    qv[qv==missing_value]=np.nan
    H[H==missing_value]=np.nan

#    print Temp
    
#    print "Length of lons {}".format(len(lons))
#    print "Length of lats {}".format(len(lats))
#    print "Length of Plevs {}".format(len(Plevs))
    
    # number of lat,lon grid nodes
    n_latitude_points = len(lats)
    n_longitude_points = len(lons)
 
    # generate the pressure levels
#    print Plevs.shape
    Pressure = np.tile(Plevs,(1,n_latitude_points,n_longitude_points,1))
    Pressure = np.transpose(Pressure,(0,2,1,3))
 
    Temp = np.transpose(Temp,(0,3,2,1))
    qv = np.transpose(qv,(0,3,2,1))
    H = np.transpose(H,(0,3,2,1))

#    print "Pressure shape {}".format(Pressure.shape)
#    print "Temp shape {}".format(Temp.shape)
#    print "qv shape {}".format(qv.shape)
#    print "H shape {}".format(H.shape)

    E = Rd/Rv
    WVapour= qv*Pressure/(E*(1-qv)+qv)
    
#    print n_longitude_points,n_latitude_points
    
    lon_vec = [i for i in range(n_longitude_points)]
    lat_vec = [i for i in range(n_latitude_points)]
    
    xx, yy = np.meshgrid(lon_vec,lat_vec)
    
#    print xx
#    print yy
#    print "xx shape {}".format(xx.shape)
#    print "yy shape {}".format(yy.shape)
    
    latgrid = np.tile(lats,(len(Plevs),n_longitude_points,1))
#    print latgrid.shape
    
    longrid = np.tile(lons,(len(Plevs),n_latitude_points,1))
#    print longrid.shape
    
    latgrid = np.transpose(latgrid,(1,2,0))
    longrid = np.transpose(longrid,(2,1,0))
    
#    print "latgrid shape {}".format(latgrid.shape)
#    print "longrid shape {}".format(longrid.shape)
    
    if sum(lons>180) > 1:
        lon0360_flag = 'y'
    else:
        lon0360_flag = 'n'
  
    Temp = np.squeeze(Temp, axis=0)
    WVapour = np.squeeze(WVapour, axis=0)
    H = np.squeeze(H, axis=0)
    Pressure = np.squeeze(Pressure, axis=0)
  
#    print "Pressure shape {}".format(Pressure.shape)
#    print "Temp shape {}".format(Temp.shape)
#    print "Wvapour shape {}".format(WVapour.shape)
#    print "H shape {}".format(H.shape)
  
    print "Temp[:,:,0] {}".format(Temp[:,:,0])
    print "Temp[:,:,1] {}".format(Temp[:,:,1])
    print "Temp[:,:,2] {}".format(Temp[:,:,2])
    print "Temp[:,:,3] {}".format(Temp[:,:,3])
    print "Temp[:,:,4] {}".format(Temp[:,:,4])
    print "Temp[:,:,5] {}".format(Temp[:,:,5])
    
    
    return (Temp,WVapour,H,Pressure,longrid,latgrid,xx,yy,lon0360_flag)

if __name__ == '__main__':

  parser = argparse.ArgumentParser(prog='aps_load_merra',
    description='Read a MERRA2 weather model file')
  parser.add_argument("wfile",help="Name of weather file to read")
#  parser.add_argument("-d","--dem",help="Name of DEM file to use (geotiff)")
#  parser.add_argument("-g","--geo_ref_file",help="Name of file for georeferencing information")
#  parser.add_argument("order_flag",help="0 - generate download file; 1 - download data")

  args = parser.parse_args()

  aps_load_merra(args.wfile)
