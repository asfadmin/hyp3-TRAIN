#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
###############################################################################
# aps_load_merra.py 
#
# Project:  APD TRAIN
# Purpose:  Load MERRA2 netcdf file
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
from netCDF4 import Dataset as NetCDFFile
import numpy as np


def aps_load_merra(wfile):
    # Set a couple of constants
    missing_value = 999999986991104
    Rv = 461.495  # [j/kg/K] water vapour
    Rd = 287.05  # [j/kg/K] specific gas constants dry air

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

    # Update no data values with NANs
    Temp[Temp == missing_value] = np.nan
    qv[qv == missing_value] = np.nan
    H[H == missing_value] = np.nan

    # number of lat,lon grid nodes
    n_latitude_points = len(lats)
    n_longitude_points = len(lons)

    # generate the pressure levels
    Pressure = np.tile(Plevs, (1, n_latitude_points, n_longitude_points, 1))
    Pressure = np.transpose(Pressure, (0, 2, 1, 3))

    Temp = np.transpose(Temp, (0, 3, 2, 1))
    qv = np.transpose(qv, (0, 3, 2, 1))
    H = np.transpose(H, (0, 3, 2, 1))
    E = Rd / Rv
    WVapour = qv * Pressure / (E * (1 - qv) + qv)

    lon_vec = [i for i in range(n_longitude_points)]
    lat_vec = [i for i in range(n_latitude_points)]
    xx, yy = np.meshgrid(lon_vec, lat_vec)

    latgrid = np.tile(lats, (len(Plevs), n_longitude_points, 1))
    longrid = np.tile(lons, (len(Plevs), n_latitude_points, 1))
    latgrid = np.transpose(latgrid, (1, 2, 0))
    longrid = np.transpose(longrid, (2, 1, 0))

    if sum(lons > 180) > 1:
        lon0360_flag = 'y'
    else:
        lon0360_flag = 'n'

    Temp = np.squeeze(Temp, axis=0)
    WVapour = np.squeeze(WVapour, axis=0)
    H = np.squeeze(H, axis=0)
    Pressure = np.squeeze(Pressure, axis=0)

    return (Temp, WVapour, H, Pressure, longrid, latgrid, xx, yy, lon0360_flag)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='aps_load_merra',
                                     description='Read a MERRA2 weather model file')
    parser.add_argument("wfile", help="Name of weather file to read")
    args = parser.parse_args()

    aps_load_merra(args.wfile)
