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
import os
import argparse
import aps_weather_model_lib as aps
from aps_merra_files import aps_merra_files
from aps_weather_model_SAR import aps_weather_model_SAR
from aps_weather_model_INSAR import aps_weather_model_INSAR
from aps_weather_model_diff import aps_weather_model_diff
import datetime
import time
import logging


def aps_weather_model(model_type, start, end, geo_ref_file=None):
    start_time = datetime.datetime.now()

    # Perform some error checking on inputs
    if start < 0 or start > 4:
        logging.error("ERROR: invalid start step {}.  Please use a number in the range of 0-4".format(start))
        exit(1)
    if end < 0 or end > 4:
        logging.error("ERROR: invalid end step {}.  Please use a number in the range of 0-4".format(end))
        exit(1)
    if end < start:
        logging.error("ERROR: End step ({}) must be greater or equal to start step ({})".format(start, end))
        exit(1)
    if geo_ref_file is not None:
        if not os.path.isfile(geo_ref_file):
            logging.error("ERROR: Unable to find geo_ref_file {}".format(geo_ref_file))
            exit(1)
    if model_type != 'era' and model_type != "merra" and model_type != "merra2":
        logging.error("ERROR: Unknown model type {}".format(model_type))
        exit(1)

    if start == 0:
        logging.info("================================================================================")
        logging.info("Step 0: Determine files to download")
        logging.info("================================================================================")
        if model_type == 'era':
            logging.error("ERROR: ERA is still under construction")
            exit(1)
        elif model_type == 'merra':
            logging.error("ERROR: MERRA is still under construction")
            exit(1)
        elif model_type == 'merra2':
            aps_merra_files(0, geo_ref_file)
        else:
            logging.error("ERROR: Unknown model type {}".format(model_type))
            exit(1)
    if start <= 1 and end >= 1:
        logging.info("================================================================================")
        logging.info("Step 1: Order and download {type} data files".format(type=model_type))
        logging.info("================================================================================")
        if model_type == 'era':
            logging.error("ERROR: ERA is still under construction")
            exit(1)
        elif model_type == 'merra':
            logging.error("ERROR: MERRA is still under construction")
            exit(1)
        elif model_type == 'merra2':
            aps_merra_files(1, geo_ref_file)
        else:
            logging.error("ERROR: Unknown model type {}".format(model_type))
            exit(1)
    if start <= 2 and end >= 2:
        logging.info("================================================================================")
        logging.info(
            "Step 2: Compute the {type} tropospheric (zenith) delay for individual dates".format(type=model_type))
        logging.info("================================================================================")
        aps_weather_model_SAR(model_type, geo_ref_file)
    if start <= 3 and end >= 3:
        logging.info("================================================================================")
        logging.info("Step 3: Compute {type} tropospheric (slant) delay for inteferograms".format(type=model_type))
        logging.info("================================================================================")
        aps_weather_model_INSAR(model_type, geo_ref_file)
    if start <= 4 and end >= 4:
        logging.info("================================================================================")
        logging.info("Step 4: Subtract delay from interferograms")
        logging.info("================================================================================")
        aps_weather_model_diff()

    logging.info("================================================================================")
    logging.info(
        "Done!!!  Total processing time {}".format(aps.timestamp(datetime.datetime.now()) - aps.timestamp(start_time)))
    logging.info("================================================================================")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='aps_weather_model',
                                     description="Calcuate tropospheric (slant) delay for a stack of interferograms: \n\
                 0 - Set up list of weather data files to download; \n\
                 1 - Download weather data; \n\
                 2 - Calculate wet and hydrostatic zenith delays; \n\
                 3 - Calculate the SAR delays; \n\
                 4 - Subtract calculated delay from interferograms")
    parser.add_argument("sstep", help="Start step to begin processing (0-4)")
    parser.add_argument("estep", help="End step to stop processing (0-4)")
    parser.add_argument("-g", "--geo_ref_file", help="Name of file for georeferencing information")
    args = parser.parse_args()

    logFile = "aps_weather_model_log.txt"
    logging.basicConfig(filename=logFile, format='%(asctime)s - %(levelname)s - %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p', level=logging.DEBUG)
    logging.getLogger().addHandler(logging.StreamHandler())
    logging.info("Starting run")

    aps_weather_model("merra2", int(args.sstep), int(args.estep), geo_ref_file=args.geo_ref_file)
