#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
###############################################################################
# aps_weather_model_SAR.py 
#
# Project:  APD TRAIN
# Purpose:  Compute SAR delays from weather model data
#           It will give the computed ZENITH hydrostatic and wet delay map in cm 
#           for the selected region.
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
from get_dem import get_ll_dem
from osgeo import gdal
import argparse
import aps_weather_model_lib as aps
import saa_func_lib as saa
import numpy as np
import scipy as sp
import scipy.integrate
import shutil
from aps_load_merra import aps_load_merra
from aps_weather_model_nan_check import aps_weather_model_nan_check
import datetime
import time
import logging
from execute import execute


# from joblib import Parallel, delayed
# import multiprocessing
# from multiprocessing import Pool


def get_DEM(geo_ref_file=None):
    region_lat_range = aps.get_range('region_lat_range', geo_ref_file=geo_ref_file)
    region_lon_range = aps.get_range('region_lon_range', geo_ref_file=geo_ref_file)

    minlat = min(region_lat_range)
    maxlat = max(region_lat_range)
    minlon = min(region_lon_range)
    maxlon = max(region_lon_range)

    if geo_ref_file is None:
        region_res = aps.get_param('region_res')
    else:
        if ".tif" in geo_ref_file:
            logging.info("Reading geotiff reference file; getting resolution.")
            (x, y, trans, proj, data) = saa.read_gdal_file(saa.open_gdal_file(geo_ref_file))
            region_res = trans[1]
        else:
            logging.error("ERROR: Unknown geo_ref_file file type {}".format(geo_ref_file))
            exit(1)

    dem_origin = aps.get_param("DEM_origin")
    dem_file = aps.get_param("DEM_file")

    if dem_origin == "asf":
        logging.info("Creating DEM file {} using the ASF DEM heap...".format(dem_file))
        get_ll_dem(minlon, minlat, maxlon, maxlat, dem_file)
    elif dem_origin == "opentopo":
        logging.info("Creating DEM file {} using opentopo...".format(dem_file))
        cmd = "wget -O%s \"http://opentopo.sdsc.edu/otr/getdem?demtype=SRTMGL1&west=%s&south=%s&east=%s&north=%s&outputFormat=GTiff\"" % (
        dem_file, minlon, minlat, maxlon, maxlat)
        execute(cmd)
    elif dem_origin == "file":
        if not os.path.isfile(dem_file):
            logging.error("ERROR: DEM file {} does not exist".format(dem_file))
            exit(1)
        else:
            logging.info("Reading in pre-existing DEM file {}".format(dem_file))
    else:
        logging.error("ERROR: Unknown DEM_origin {} read from file parms_aps.txt".format(dem_origin))
        exit(1)

    (x, y, trans, proj, data) = saa.read_gdal_file(saa.open_gdal_file(dem_file))
    if np.abs(float(trans[1]) - float(region_res)) > 10e-9:
        logging.info("Resampling DEM file to match region_res")
        gdal.Warp("tmp.tif", dem_file, xRes=region_res, yRes=region_res, resampleAlg="cubic", dstNodata=-32767,
                  format="GTiff")
        shutil.move("tmp.tif", dem_file)
        (x, y, trans, proj, data) = saa.read_gdal_file(saa.open_gdal_file(dem_file))
    nncols = x
    nnrows = y
    demres = trans[1]
    ullat = trans[3]
    lrlat = trans[3] + y * trans[5]
    ullon = trans[0]
    lrlon = trans[0] + x * trans[1]
    xmin = min(ullon, lrlon)
    xmax = max(ullon, lrlon)
    ymin = min(ullat, lrlat)
    ymax = max(ullat, lrlat)

    return (data, minlon, maxlon, minlat, maxlat, demres, nncols, nnrows)


def write_out_file(name, xi, yi, data):
    fid = open(name, 'wb')
    data_write = np.array(
        [np.reshape(xi, -1, order='F'), np.reshape(yi, -1, order='F'), np.reshape(data, -1, order='F')])
    data_write = np.reshape(data_write, (3, -1))
    data_write = np.transpose(data_write)
    data_write.tofile(fid)
    fid.close()


def inter2d(xivec, yivec, lonlist_matrix, latlist_matrix, cdstack, n):
    tck = sp.interpolate.bisplrep(lonlist_matrix, latlist_matrix, cdstack[:, :, n], s=3)
    newslice = sp.interpolate.bisplev(yivec, xivec, tck)
    return newslice


def aps_weather_model_SAR(model_type, geo_ref_file=None):
    start = datetime.datetime.now()

    # Constants
    # Parameters from Jolivet et al 2011, GRL
    Rd = 287.05  # [j/kg/K] specific gas constants dry air
    Rv = 461.495  # [j/kg/K] water vapour
    k1 = 0.776  # [K/Pa]
    k2 = 0.716  # [K/Pa]
    k3 = 3.75e3  # [K^2/Pa]
    g0 = 9.80665
    Rmax = 6378137
    Rmin = 6356752

    zref = 15000  # zref for integration- tropopause
    zincr = 15  # z increment spacing for vertical interpolation before integration
    vertres = 100  # vertical resolution for comparing dem to vert profiles of delay

    XI = [int(i) for i in range(0, zref + 1, zincr)]

    path = aps.get_param('{}_datapath'.format(model_type))
    utc = float(aps.get_param('UTC_sat'))
    dates, tmp = aps.get_date_list()
    n_dates = len(dates)

    if geo_ref_file is not None:
        if not os.path.isfile(geo_ref_file):
            logging.error("ERROR: Unable to find geo_ref_file {}".format(geo_ref_file))
            exit(1)

    dem, xmin, xmax, ymin, ymax, smpres, nncols, nnrows = get_DEM(geo_ref_file)

    lonmin = np.floor(xmin) - 1
    lonmax = np.ceil(xmax) + 1
    latmin = np.floor(ymin) - 1
    latmax = np.ceil(ymax) + 1

    maxdem = np.ceil(np.max(np.max(dem)) / 100) * 100 + 100
    cdslices = int(maxdem / vertres) + 1
    cdI = [x for x in range(0, int(maxdem) + 1, vertres)]

    rounddem = np.round(dem / vertres)
    rounddem[dem < 0] = 0
    rounddem[np.isnan(dem)] = 0
    del dem

    if cdslices != len(cdI):
        logging.error("INTERNAL ERROR: length mistmatch")
        logging.error("cdslices {}".format(cdslices))
        logging.error("len(cdI) {}".format(len(cdI)))
        logging.error("cdI {}".format(cdI))
        exit(1)

    logging.info("Interpolate to a maximum dem height of {HEIGHT}".format(HEIGHT=maxdem))
    logging.info("Using {} slices".format(cdslices))

    date, time, frac = aps.times(utc, dates)

    input_file_names = aps.file_names(model_type, date, time, path)
    length = len(input_file_names) / 2

    for d in range(length):
        lasttime = datetime.datetime.now()
        no_data = 0
        for kk in range(2):
            if kk == 0:
                wfile = input_file_names[d]
            else:
                wfile = input_file_names[d + length]

            if not os.path.isfile(wfile):
                no_data = no_data + 1
            else:
                logging.info("Loading input weather file: {}".format(wfile))

                # Load the weather model data
                if model_type == "merra2":
                    (Temp, e, H, P, longrid, latgrid, xx, yy, lon0360_flag) = aps_load_merra(wfile)
                else:
                    logging.info("Model type {} not implemented!".format(model_type))
                    exit(1)

                # deal with NANs
                (Temp, e) = aps_weather_model_nan_check(Temp, e, P, longrid, latgrid)

                # define weather model grid nodes
                latlist = np.reshape(latgrid[:, :, 1], -1)
                lonlist = np.reshape(longrid[:, :, 1], -1)
                xlist = np.reshape(xx, -1, order='F')
                ylist = np.reshape(yy, -1, order='F')
                lat_res = np.abs(np.diff(np.unique(latgrid))) * 1.5
                lat_res = lat_res[0]
                lon_res = np.abs(np.diff(np.unique(longrid))) * 1.5
                lon_res = lon_res[0]

                if lon0360_flag == 'y':
                    if xmin < 0:
                        xmin = xmin + 360
                    if xmax < 0:
                        xmax = xmax + 360

                # ix = np.where(ymin-lat_res<= latlist and latlist <= ymax+lat_res and xmin-lon_res<= lonlist and lonlist<= xmax+lon_res)
                ix1 = np.where(ymin - lat_res <= latlist)
                ix2 = np.where(latlist <= ymax + lat_res)
                ix3 = np.intersect1d(ix1, ix2)
                ix1 = np.where(xmin - lon_res <= lonlist)
                ix2 = np.where(lonlist <= xmax + lon_res)
                ix4 = np.intersect1d(ix1, ix2)
                ix = np.intersect1d(ix3, ix4)

                xlist = xlist[ix]
                ylist = ylist[ix]
                latlist = latlist[ix]
                lonlist = lonlist[ix]

                ulatlist = np.unique(latlist)
                ulonlist = np.unique(lonlist)
                numy = int(len(ulatlist))
                numx = int(len(ulonlist))
                uxlist = np.unique(xlist)
                uylist = np.unique(ylist)

                # map of g with latitude
                g = 9.80616 * (1. - 0.002637 * np.cos(2 * np.deg2rad(latgrid)) + 0.0000059 * np.square(
                    np.cos(2. * np.deg2rad(latgrid))))

                # map of Re with latitude
                Re = np.sqrt(1. / ((np.square(np.cos(np.deg2rad(latgrid))) / np.square(Rmax)) + (
                            np.square(np.sin(np.deg2rad(latgrid))) / np.square(Rmin))))

                # Calculate Geometric Height, Z
                Z = (H * Re) / (g / g0 * Re - H)

                midx = int(round(np.mean(uxlist)))
                midy = int(round(np.mean(uylist)))

                glocal = g[midx, midy, 0]
                Rlocal = Re[midx, midy, 0]

                del g
                del Re
                del H

                cdstack_dry = np.zeros((numy, numx, cdslices), dtype=np.float32)
                cdstack_wet = np.zeros((numy, numx, cdslices), dtype=np.float32)

                # Interpolate Temp P and e from 0:20:15000 m
                # then integrate using trapz to estimate delay as function of height
                for i in range(numx):
                    for j in range(numy):
                        xn = uxlist[i]
                        yn = uylist[j]

                        # interpolate at res zincr before doing integration
                        X = np.squeeze(Z[xn, yn, :])
                        Ye = np.squeeze(e[xn, yn, :])
                        f = sp.interpolate.splrep(X, Ye)
                        YeI = sp.interpolate.splev(XI, f, der=0)
                        YeI = YeI * 100

                        Yp = np.squeeze(P[xn, yn, :])
                        f = sp.interpolate.splrep(X, Yp)
                        YPI = sp.interpolate.splev(XI, f, der=0)
                        YPI = YPI * 100

                        YT = np.squeeze(Temp[xn, yn, :])
                        f = sp.interpolate.splrep(X, YT)
                        YTI = sp.interpolate.splev(XI, f, der=0)

                        tmp1 = (k2 - (Rd * k1 / Rv)) * YeI / YTI + k3 * YeI / (YTI * YTI)
                        Lw = (10 ** -6) * -1 * np.flipud(
                            sp.integrate.cumtrapz(np.flipud(tmp1), np.flipud(XI), initial=0))

                        f = sp.interpolate.splrep(XI, Lw)
                        LwI = sp.interpolate.splev(cdI, f, der=0)
                        Ld = 10 ** -6 * ((k1 * Rd / glocal) * (YPI - YPI[zref / zincr]))
                        f = sp.interpolate.splrep(XI, Ld)
                        LdI = sp.interpolate.splev(cdI, f, der=0)

                        cdstack_dry[j, i, :] = LdI
                        cdstack_wet[j, i, :] = LwI

                del LwI, LdI

                xsmpres = (xmax - xmin) / nncols
                ysmpres = (ymax - ymin) / nnrows

                xivec = np.linspace(xmin + 0.5 * xsmpres, xmax - 0.5 * xsmpres, nncols)
                yivec = np.linspace(ymin + 0.5 * ysmpres, ymax - 0.5 * ysmpres, nnrows)
                [xi, yi] = np.meshgrid(xivec, yivec)

                ix_temp = np.diff(lonlist)
                ix_temp = np.where(ix_temp != 0)
                ixi_temp = ix_temp[0] + 1
                lonlist_matrix = np.reshape(lonlist, (ixi_temp[0], -1), order='F')
                latlist_matrix = np.reshape(latlist, (ixi_temp[0], -1), order='F')
                del ix_temp
                del ixi_temp

                cdstack_interp_dry = np.zeros((nnrows, nncols, cdslices), dtype=np.float32)
                logging.info("processing dry stack")
                for n in range(cdslices):
                    #                    f = scipy.interpolate.interp2d(lonlist_matrix,latlist_matrix,cdstack_dry[:,:,n])
                    #                    newslice = f(xivec,yivec)

                    tck = sp.interpolate.bisplrep(lonlist_matrix, latlist_matrix, cdstack_dry[:, :, n], s=3)
                    newslice = sp.interpolate.bisplev(yivec, xivec, tck)

                    cdstack_interp_dry[:, :, n] = np.flipud(newslice)
                del cdstack_dry

                #                Parallel(n_jobs=8)(delayed(inter2d)(xivec,yivec,lonlist_matrix,latlist_matrix,cdstack,cdstack_interp_dry,n)
                #                           for n in range(cdslices))

                #                pool = Pool(processes=8)
                #                cdstack_interp_dry = [pool.apply(inter2d,args=(xivec,yivec,lonlist_matrix,latlist_matrix,cdstack,n)) for n in range(cdslices)]

                cdstack_interp_wet = np.zeros((nnrows, nncols, cdslices), dtype=np.float32)
                logging.info("processing wet stack")

                lonlist_matrix = lonlist_matrix[0, :]
                latlist_matrix = latlist_matrix[:, 0]

                for n in range(cdslices):
                    #                    f = sp.interpolate.interp2d(lonlist_matrix,latlist_matrix,cdstack_wet[:,:,n])

                    f = sp.interpolate.RectBivariateSpline(latlist_matrix, lonlist_matrix, cdstack_wet[:, :, n], kx=1,
                                                           ky=1, s=0)
                    newslice = f(yivec, xivec)

                    #                    tck = sp.interpolate.bisplrep(lonlist_matrix,latlist_matrix,cdstack_wet[:,:,n],s=5)
                    #                    newslice = sp.interpolate.bisplev(yivec,xivec,tck)

                    cdstack_interp_wet[:, :, n] = np.flipud(newslice)
                del cdstack_wet

                # keeping the coordinates in the same grid as the data
                xi = np.flipud(xi)
                yi = np.flipud(yi)

                # Pull out delays from cdstack layers that match dem heights
                wetcorrection = np.ones((nnrows, nncols), dtype=np.float32)
                hydrcorrection = np.ones((nnrows, nncols), dtype=np.float32)

                for i in range(nnrows):
                    for j in range(nncols):
                        wetcorrection[i, j] = cdstack_interp_wet[i, j, int(rounddem[i, j])]
                        hydrcorrection[i, j] = cdstack_interp_dry[i, j, int(rounddem[i, j])]

                if kk == 0:
                    wetcorrection1 = wetcorrection
                    drycorrection1 = hydrcorrection
                elif kk == 1:
                    wetcorrection2 = wetcorrection
                    drycorrection2 = hydrcorrection

                del wetcorrection
                del hydrcorrection
                del cdstack_interp_wet
                del cdstack_interp_dry

        if no_data == 0:

            logging.info("Writing output files")

            outfile_wet_before = path + "/" + date[d] + "/" + date[d] + '_ZWD_before.xyz'
            outfile_wet_after = path + "/" + date[d] + "/" + date[d] + '_ZWD_after.xyz'

            outfile_dry_before = path + "/" + date[d] + "/" + date[d] + '_ZHD_before.xyz'
            outfile_dry_after = path + "/" + date[d] + "/" + date[d] + '_ZHD_after.xyz'

            write_out_file(outfile_wet_before, xi, yi, wetcorrection1)
            write_out_file(outfile_wet_after, xi, yi, wetcorrection2)

            write_out_file(outfile_dry_before, xi, yi, drycorrection1)
            write_out_file(outfile_dry_after, xi, yi, drycorrection2)

            # Output wet correction
            wetcorrection = wetcorrection1 * frac[d] + wetcorrection2 * frac[d + length]
            del wetcorrection1
            del wetcorrection2
            wetcorrection = wetcorrection * 100  # delay in CM
            outfile = path + "/" + date[d] + "/" + date[d] + '_ZWD.xyz'
            write_out_file(outfile, xi, yi, wetcorrection)
            del wetcorrection

            # Output hydrostatic correction
            hydrcorrection = drycorrection1 * frac[d] + drycorrection2 * frac[d + length]
            del drycorrection1
            del drycorrection2
            hydrcorrection = hydrcorrection * 100  # delay in CM
            hydroutfile = path + "/" + date[d] + "/" + date[d] + '_ZHD.xyz'
            write_out_file(hydroutfile, xi, yi, hydrcorrection)
            del hydrcorrection

            elapsed = aps.timestamp(datetime.datetime.now()) - aps.timestamp(lasttime)
            logging.info("{D} completed out of {L} in time {S}".format(D=d + 1, L=length, S=elapsed))
        else:
            logging.info("{D} completed out of {L} (NO DATA)".format(D=d + 1, L=length))

    elapsed = aps.timestamp(datetime.datetime.now()) - aps.timestamp(start)
    logging.info("Done!!! Completed all files in time {}".format(elapsed))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='aps_weather_model_SAR',
                                     description='Calculate zenith wet and hydrostatic delays')
    parser.add_argument("-g", "--geo_ref_file", help="Name of file for georeferencing information")
    args = parser.parse_args()

    logFile = "aps_weather_model_SAR_log.txt"
    logging.basicConfig(filename=logFile, format='%(asctime)s - %(levelname)s - %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p', level=logging.DEBUG)
    logging.getLogger().addHandler(logging.StreamHandler())
    logging.info("Starting run")

    aps_weather_model_SAR("merra2", geo_ref_file=args.geo_ref_file)
