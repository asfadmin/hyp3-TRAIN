#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
###############################################################################
# aps_weather_model_nan_check.py 
#
# Project:  APD TRAIN
# Purpose:  
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
import numpy as np
import logging
from scipy.interpolate import griddata

def aps_weather_model_nan_check(Temp,e,P,longrid,latgrid):

    Pressure_level = np.squeeze(P[1,1,:])
    lon = longrid[:,:,1]
    lat = latgrid[:,:,1]
    step_level = [i for i in range(len(Pressure_level)-1,-1,-1)]
    if Pressure_level[0]<Pressure_level[-1]:
        step_level = [i for i in range(len(Pressure_level))]
    dims = Temp.shape
    n_pixels = dims[0] * dims[1]
    grid = np.zeros((n_pixels,2))
    grid[:,0] = lon.flatten()
    grid[:,1] = lat.flatten()        

    for k in range(len(Pressure_level)):
        test = np.squeeze(Temp[:,:,step_level[k]])
   	ix_nan = np.isnan(test)
        count = np.zeros((dims[0],dims[1]))
        for i in range(dims[0]):
           for j in range(dims[1]):
               if ix_nan[i,j] == False:
                   count[i,j] += int(1)
               else:
                    ix_nan[i,j] = True
        n_nan = n_pixels - sum(sum(count))        
        if n_nan < n_pixels-3:
            points = np.zeros((len(lon[~ix_nan]),2))
            points[:,0] = np.reshape(lon[~ix_nan],-1,1)
            points[:,1] = np.reshape(lat[~ix_nan],-1,1)
            temp_data = Temp[:,:,step_level[k]]
            data = np.reshape(temp_data[~ix_nan],-1,1)
            temp_data = griddata(points,data,grid,'nearest')
            Temp[:,:,step_level[k]] = np.reshape(temp_data,(dims[0],dims[1]))
            temp_data = e[:,:,step_level[k]]
            data = np.reshape(temp_data[~ix_nan],-1,1)
            temp_data = griddata(points,data,grid,'nearest')
            e[:,:,step_level[k]] = np.reshape(temp_data,(dims[0],dims[1]))
        elif n_nan == n_pixels:
            if k==1:
                logging.warning("WARNING:  The top of the atmosphere has no data")
            Temp[:,:,step_level[k]]=Temp[:,:,step_level[k-1]]
            e[:,:,step_level[k]]=e[:,:,step_level[k-1]]
        else:
            Temp[:,:,step_level[k]]=np.nanmean(Temp[:,:,step_level[k]])
            e[:,:,step_level[k]]=np.nanmean(e[:,:,step_level[k]])
    return  (Temp,e)
    

