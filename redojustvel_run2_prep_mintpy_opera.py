#!/usr/bin/env python
"""Truncate in time existing OPERA DISP-S1 derived time-series.

This script allows users to:
1. Truncate existing time-series with specified temporal subsets
2. Regenerate velocity maps with the time-series subset

Example:
    redojustvel_run2_prep_mintpy_opera.py -u "outputs/*.nc" -m static_lyrs
        --geom-dir geometry -o mintpy_output -s 20160701 -e 20200901

Input Requirements:
    - DISP-S1 NetCDF files containing displacements (needed to obtain metadata)
    - Static layer files with geometry information
    - Specify temporal subsetting with a start and/or end-date

Outputs:
    - timeseries.h5: Cumulative displacement time series
    - velocity.h5: Linear displacement rates

Dependencies:
    mintpy, gdal, rasterio, h5py, numpy, pandas, cartopy
    opera_utils (for handling OPERA-specific file formats)

Base code
Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi 
Author: Talib Oliver Cabrerra, Scott Staniewicz 

DISP-S1 Implementation
Author: Simran S Sangha, Jinwoo Kim
February, 2025
"""

# Standard library imports
import argparse
import glob
import os
import sys
import warnings
from datetime import datetime as dt
from pathlib import Path
from typing import Sequence
import time

# Third-party imports
import asf_search as asf
import h5py
import networkx as nx
import numpy as np
import pandas as pd
import rasterio
from osgeo import gdal
from packaging.version import Version
from tqdm import tqdm
import psutil

# Suppress warnings
warnings.filterwarnings("ignore")

# Add the src directory to sys.path
sys.path.append(str(Path(__file__).parent / "src"))

# Local application/library-specific imports
from mintpy.cli import (
    generate_mask,
    mask,
)
from mintpy.utils import (
    arg_utils,
    ptime,
    readfile,
    writefile,
)
from mintpy.utils.utils0 import (
    azimuth2heading_angle,
    calc_azimuth_from_east_north_obs,
)
from pst_dolphin_utils import (
    BackgroundRasterWriter,
    datetime_to_float,
    estimate_velocity,
    full_suffix,
    get_raster_crs,
    get_raster_gt,
    get_raster_xysize,
    HDF5StackReader,
    process_blocks,
)

OPERA_DATASET_ROOT = './'

EXAMPLE = """example:
    redojustvel_run2_prep_mintpy_opera.py -u "outputs/*.nc" -m static_lyrs
        --geom-dir geometry -o mintpy_output -s 20160701 -e 20200901
"""

def _create_parser():
    parser = argparse.ArgumentParser(
        description="Prepare Sweets products for MintPy",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=EXAMPLE,
    )

    parser.add_argument(
        "-u",
        "--unw-file-glob",
        type=str,
        default="./interferograms/unwrapped/*.unw.tif",
        help="path pattern of unwrapped interferograms (default: %(default)s).",
    )
    parser.add_argument(
        "-g",
        "--geom-dir",
        default="./geometry",
        help="Geometry directory (default: %(default)s).",
    )
    parser.add_argument(
        "-m",
        "--meta-file",
        type=str,
        help="GSLC metadata file or directory",
    )
    parser.add_argument(
        "-s",
        "--start-date",
        dest='startDate',
        default=None,
        help="remove/drop interferograms with date earlier than "
             "start-date in YYMMDD or YYYYMMDD format",
    )
    parser.add_argument(
        "-e",
        "--end-date",
        dest='endDate',
        default=None,
        help="remove/drop interferograms with date later than "
             "end-date in YYMMDD or YYYYMMDD format",
    )
    parser.add_argument(
        "-o",
        "--out-dir",
        type=str,
        default="./mintpy",
        help="output directory (default: %(default)s).",
    )
    parser.add_argument(
        "-r",
        "--range",
        dest="lks_x",
        type=int,
        default=1,
        help=(
            "number of looks in range direction, for multilooking applied after fringe"
            " processing.\nOnly impacts metadata. (default: %(default)s)."
        ),
    )
    parser.add_argument(
        "-a",
        "--azimuth",
        dest="lks_y",
        type=int,
        default=1,
        help=(
            "number of looks in azimuth direction, for multilooking applied after"
            " fringe processing.\nOnly impacts metadata. (default: %(default)s)."
        ),
    )
    parser.add_argument(
        "--zero-mask",
        dest="zero_mask",
        action="store_true",
        help="Mask all pixels with zero value in unw phase",
    )
    parser.add_argument(
        "--work-dir",
        dest="work_dir",
        type=str,
        default="./raider_intermediate_workdir",
        help="Working directory for tropospheric correction intermediates"
    )

    parser = arg_utils.add_subset_argument(parser, geo=True)

    return parser


def cmd_line_parse(iargs=None):
    """Create the command line parser"""
    parser = _create_parser()
    inps = parser.parse_args(args=iargs)

    # in case meta_file is input as wildcard
    inps.meta_file = sorted(glob.glob(inps.meta_file))[0]

    return inps


def prepare_metadata(meta_file, int_file, geom_dir, nlks_x=1, nlks_y=1):
    """Get the metadata from the GSLC metadata file and the unwrapped interferogram."""
    print("-" * 50)

    cols, rows = get_raster_xysize(int_file)

    meta_compass = h5py.File(meta_file, "r")
    meta = {}

    geotransform = get_raster_gt(int_file)
    meta["LENGTH"] = rows
    meta["WIDTH"] = cols

    meta["X_FIRST"] = geotransform[0]
    meta["Y_FIRST"] = geotransform[3]
    meta["X_STEP"] = geotransform[1]
    meta["Y_STEP"] = geotransform[5]
    meta["X_UNIT"] = meta["Y_UNIT"] = "meters"

    crs = get_raster_crs(int_file)
    meta["EPSG"] = crs.to_epsg()

    if str(meta["EPSG"]).startswith('326'):
         meta["UTM_ZONE"] = str(meta["EPSG"])[3:] + 'N'
    else:
         meta["UTM_ZONE"] = str(meta["EPSG"])[3:] + 'S'

    if "/science" in meta_compass:
        root = "/science/SENTINEL1/CSLC"
        processing_ds = f"{root}/metadata/processing_information"
        burst_ds = f"{processing_ds}/s1_burst_metadata"
        if burst_ds not in meta_compass:
            burst_ds = f"{processing_ds}/input_burst_metadata"
    else:
        root = OPERA_DATASET_ROOT
        processing_ds = f"{root}/metadata/processing_information"
        burst_ds = f"{processing_ds}/input_burst_metadata"

    meta["WAVELENGTH"] = meta_compass[f"{burst_ds}/wavelength"][()]
    meta["RANGE_PIXEL_SIZE"] = meta_compass[f"{burst_ds}/range_pixel_spacing"][()]
    meta["AZIMUTH_PIXEL_SIZE"] = 14.1
    meta["EARTH_RADIUS"] = 6371000.0

    # get heading from azimuth angle
    geom_path = Path(geom_dir)
    file_to_path = {
        "los_east": geom_path / "los_east.tif",
        "los_north": geom_path / "los_north.tif",
    }
    dsDict = {}
    for dsName, fname in file_to_path.items():
        data = readfile.read(fname, datasetName=dsName)[0]
        data[data == 0] = np.nan
        dsDict[dsName] = data
    azimuth_angle, _, _ = get_azimuth_ang(dsDict)
    azimuth_angle = np.nanmean(azimuth_angle)
    heading = azimuth2heading_angle(azimuth_angle)
    meta["HEADING"] = heading

    t0 = dt.strptime(
        meta_compass[f"{burst_ds}/sensing_start"][()].decode("utf-8"),
        "%Y-%m-%d %H:%M:%S.%f",
    )
    t1 = dt.strptime(
        meta_compass[f"{burst_ds}/sensing_stop"][()].decode("utf-8"),
        "%Y-%m-%d %H:%M:%S.%f",
    )
    t_mid = t0 + (t1 - t0) / 2.0
    meta["CENTER_LINE_UTC"] = (
        t_mid - dt(t_mid.year, t_mid.month, t_mid.day)
    ).total_seconds()
    meta["HEIGHT"] = 750000.0
    meta["STARTING_RANGE"] = meta_compass[f"{burst_ds}/starting_range"][()]
    meta["PLATFORM"] = meta_compass[f"{burst_ds}/platform_id"][()].decode("utf-8")
    meta["ORBIT_DIRECTION"] = meta_compass[f"{root}/metadata/orbit/orbit_direction"][
        ()
    ].decode("utf-8")
    meta["ALOOKS"] = 1
    meta["RLOOKS"] = 1

    # apply optional user multilooking
    if nlks_x > 1:
        meta["RANGE_PIXEL_SIZE"] = str(float(meta["RANGE_PIXEL_SIZE"]) * nlks_x)
        meta["RLOOKS"] = str(float(meta["RLOOKS"]) * nlks_x)

    if nlks_y > 1:
        meta["AZIMUTH_PIXEL_SIZE"] = str(float(meta["AZIMUTH_PIXEL_SIZE"]) * nlks_y)
        meta["ALOOKS"] = str(float(meta["ALOOKS"]) * nlks_y)

    return meta


def _get_xy_arrays(atr):
    x0 = float(atr["X_FIRST"])
    y0 = float(atr["Y_FIRST"])
    x_step = float(atr["X_STEP"])
    y_step = float(atr["Y_STEP"])
    rows = int(atr["LENGTH"])
    cols = int(atr["WIDTH"])
    x_arr = x0 + x_step * np.arange(cols)
    y_arr = y0 + y_step * np.arange(rows)
    # Shift by half pixel to get the centers
    x_arr += x_step / 2
    y_arr += y_step / 2

    return x_arr, y_arr


def _get_date_pairs(filenames):
    str_list = [Path(f).stem for f in filenames]
    basenames_noext = [str(f).replace(full_suffix(f), "") for f in str_list]

    date_pairs = []
    for i in basenames_noext:
        num_parts = i.split('_')
        if len(num_parts) == 9:
            date_pair = f'{num_parts[6][:8]}_{num_parts[7][:8]}'
            date_pairs.append(date_pair)

    return date_pairs


def get_azimuth_ang(dsDict):
    """Compute the azimuth angle from east/north coefficients"""
    east = dsDict["los_east"]
    north = dsDict["los_north"]
    azimuth_angle = calc_azimuth_from_east_north_obs(east, north)

    return azimuth_angle, east, north


def main(iargs=None):
    """Run the preparation functions."""
    inps = cmd_line_parse(iargs)

    start_time = time.time()

    product_files = sorted(
        glob.glob(inps.unw_file_glob),
        key=lambda x: dt.strptime(
            x.split('_')[-3][:8], '%Y%m%d'
        )
    )

    # track product version
    track_version = []
    for i in product_files:
        fname = os.path.basename(i)
        version_n = Version(fname.split('_')[-2].split('v')[1])
        track_version.append(version_n)

    # exit if multiple versions are found
    track_version = list(set(track_version))
    if len(track_version) > 1:
        raise Exception(f'Multiple file version increments ({track_version}) '
                        'found in specified input. Version increments are ' 
                        'not compatible. '
                        'delete the PDF.')

    # pass unw conversion factor, which depends on the product version
    track_version = track_version[0]
    if track_version == Version('0.3'):
        disp_lyr_name = 'unwrapped_phase'
    if track_version >= Version('0.4'):
        disp_lyr_name = 'displacement'

    # append appropriate NETCDF prefixes
    unw_files = \
        [f'NETCDF:"{i}":{disp_lyr_name}' for i in product_files]

    # get desired start/end time(s), if specified
    if inps.startDate is not None:
        startDate = int(inps.startDate)

    if inps.endDate is not None:
        endDate = int(inps.endDate)

    # check static layer naming convention
    static_dir = Path(inps.meta_file)
    allcaps_geometry = True
    static_files = sorted(Path(static_dir).glob("*STATIC_*.h5"))
    # capture alternate filename convention
    if static_files == []:
        allcaps_geometry = False

    # translate input options
    # metadata
    meta_file = Path(inps.meta_file)
    if meta_file.is_dir():
        # Search for the line of sight static_layers file
        try:
            # Grab the first one in in the directory
            if allcaps_geometry is True:
                meta_file = next(meta_file.rglob("*STATIC_*.h5"))
            else:
                meta_file = next(meta_file.rglob("static_*.h5"))
        except StopIteration:
            raise ValueError(f"No static layers file found in {meta_file}")

    meta = prepare_metadata(
        meta_file, unw_files[0], geom_dir=inps.geom_dir,
        nlks_x=inps.lks_x, nlks_y=inps.lks_y
    )

    # prepare TS file
    og_ts_file = os.path.join(inps.out_dir, "timeseries.h5")
    geom_file = os.path.join(inps.out_dir, "geometryGeo.h5")
    all_outputs = [og_ts_file]
    ref_meta = readfile.read_attribute(og_ts_file)

    # generate velocity fit(s)
    ts_dict = {}
    ts_dict['velocity'] = og_ts_file

    shortwvl_lyrs_path = os.path.join(inps.out_dir,
            'short_wavelength_displacement.h5')
    if os.path.exists(shortwvl_lyrs_path):
        ts_dict['velocity_shortwvl'] = shortwvl_lyrs_path

    dem_error_path = os.path.join(inps.out_dir,
        "timeseries_demErr.h5")
    if os.path.exists(dem_error_path):
        ts_dict['velocity_demErr'] = dem_error_path

    era5_corr_path = os.path.join(inps.out_dir,
        "timeseries_ERA5.h5")
    if os.path.exists(era5_corr_path):
        ts_dict['velocity_ERA5'] = era5_corr_path

    era5_demErr_corr_path = os.path.join(inps.out_dir,
        "timeseries_ERA5_demErr.h5")
    if os.path.exists(era5_demErr_corr_path):
        ts_dict['velocity_ERA5_demErr'] = era5_demErr_corr_path

    # check if files need to be truncated
    if inps.startDate is not None or inps.endDate is not None:
        for ts_name in ts_dict.values():
            # get list of times WRT to the reference time
            # and also pass the TS data
            with h5py.File(ts_name, 'r') as f:
                ts_date_list = f['date'][:]
                og_ts_date_list = [int(i) for i in ts_date_list]
                ts_date_list = [int(i) for i in ts_date_list]

            # filter dates
            if inps.startDate is not None:
                ts_date_list = [i for i in ts_date_list if i > startDate]
                ts_date_list.sort()
            if inps.endDate is not None:
                ts_date_list = [i for i in ts_date_list if i < endDate]
                ts_date_list.sort()

            # Only proceed if specified temporal sampling would actually
            # truncate the TS
            if len(og_ts_date_list) == len(ts_date_list):
                print(
                    f"Specified start date {inps.startDate} "
                    f"and end date {inps.endDate} do not "
                    f"change, no modification to {ts_name} necessary"
                )
                continue

            # rename orig TS file
            ts_name_base = os.path.basename(ts_name)
            full_ts_name = os.path.join(inps.out_dir, f'full_{ts_name_base}')
            os.rename(ts_name, full_ts_name)

            # define dataset structure
            dates = np.array(ts_date_list, dtype=np.string_)
            num_date = len(dates)
            pbase = np.zeros(num_date, dtype=np.float32)
            row = int(ref_meta['LENGTH'])
            col = int(ref_meta['WIDTH'])
            ds_name_dict = {
                "date": [dates.dtype, (num_date,), dates],
                "bperp": [np.float32, (num_date,), pbase],
                "timeseries": [np.float32, (num_date, row, col), None],
            }

            # Initialize HDF5 file
            writefile.layout_hdf5(ts_name, ds_name_dict, metadata=ref_meta)

            # Initialize with zeros
            with h5py.File(ts_name, "r+") as f:
                f["timeseries"][0] = np.zeros((row, col), dtype=np.float32)
            
            # get reference data slice
            # to subtract from the other dates
            ref_date = ts_date_list[0]
            ref_data_slice, _ = readfile.read(
                full_ts_name, datasetName=f'timeseries-{ref_date}')

            # loop through all other layers
            for ts_ind, ts_date in enumerate(ts_date_list[1:]):
                data_slice, _ = readfile.read(
                    full_ts_name, datasetName=f'timeseries-{ts_date}')
                # Write data slice to file
                with h5py.File(ts_name, "r+") as f:
                    f["timeseries"][ts_ind+1] = data_slice - ref_data_slice

    for vel_name, ts_name in ts_dict.items():
        # first set variables
        dolphin_ref_tif = os.path.join(inps.out_dir, 'dolphin_reference.tif')
        dolphin_vel_file = os.path.join(inps.out_dir, f'{vel_name}.tif')
        vel_file = os.path.join(inps.out_dir, f'{vel_name}.h5')
        keep_open = False
        dset_names = 'timeseries'
        num_threads = 6
        block_shape = (256, 256)

        # rename orig vel files
        for vel_output in [dolphin_vel_file, vel_file]:
            if os.path.exists(vel_output):
                vel_name_base = os.path.basename(vel_output)
                full_vel_name = os.path.join(inps.out_dir,
                    f'full_{vel_name_base}')
                os.rename(vel_output, full_vel_name)

        # extract one product to serve as a reference file
        ds = gdal.Translate(dolphin_ref_tif, unw_files[0])
        ds = None

        # get list of times WRT to the reference time
        # and also pass the TS data
        with h5py.File(ts_name, 'r') as f:
            ts_date_list = f['date'][:]
            ts_data = f['timeseries'][:]

        x_arr = [
            dt.strptime(
                date.decode('utf-8'), '%Y%m%d'
            )
            for date in ts_date_list
        ]
        x_arr = datetime_to_float(x_arr)

        # initiate dolphin file object
        writer = BackgroundRasterWriter(dolphin_vel_file,
            like_filename=dolphin_ref_tif)

        # run dolphin velocity fitting algorithm in blocks
        def read_and_fit(
            readers: Sequence[HDF5StackReader],
            rows: slice, cols: slice
        ) -> tuple[np.ndarray, slice, slice]:

            # Only use the cor_reader if it's the same shape as the unw_reader
            if len(readers) == 2:
                unw_reader, cor_reader = readers
                unw_stack = unw_reader[:, rows, cols]
                weights = cor_reader[:, rows, cols]
                cor_threshold = 0.4
                weights[weights < cor_threshold] = 0
            else:
                unw_stack = readers[0][:, rows, cols]
                weights = None

            # Fit a line to each pixel with weighted least squares
            return (
                estimate_velocity(
                    x_arr=x_arr,
                    unw_stack=unw_stack,
                    weight_stack=weights
                ),
                rows,
                cols,
            )

        readers = [ts_data]
        process_blocks(
            readers=readers,
            writer=writer,
            func=read_and_fit,
            block_shape=block_shape,
            num_threads=num_threads,
        )

        writer.notify_finished()

        # delete temporary reference file
        os.remove(dolphin_ref_tif)

        # update metadata field
        with h5py.File(ts_name, 'r') as f:
            ts_date_list = f['date'][:]
            start_date = str(int(ts_date_list[0]))
            end_date = str(int(ts_date_list[-1]))

        meta["DATA_TYPE"] = 'float32'
        meta["DATE12"] =  start_date + '_' + end_date
        meta["FILE_PATH"] = ts_name
        meta["FILE_TYPE"] = 'velocity'
        meta["NO_DATA_VALUE"] = 'none'
        meta["PROCESSOR"] = 'dolphin'
        meta["REF_DATE"] = start_date
        meta["START_DATE"] = start_date
        meta["END_DATE"] = end_date
        meta["UNIT"] = 'm/year'
        # apply reference point to velocity file
        if ref_meta is not None:
            meta['REF_LAT'] = ref_meta['REF_LAT']
            meta['REF_LON'] = ref_meta['REF_LON']
            meta['REF_Y'] = ref_meta['REF_Y']
            meta['REF_X'] = ref_meta['REF_X']

        # initiate HDF5 file
        row = int(meta['LENGTH'])
        col = int(meta['WIDTH'])
        ds_name_dict = {
            "velocity": [np.float32, (row, col), None],
        }
        writefile.layout_hdf5(vel_file, ds_name_dict, metadata=meta)

        # writing data to HDF5 file
        print("writing data to HDF5 file {} with a mode ...".format(vel_file))
        with h5py.File(vel_file, "a") as f:
            vel_arr = gdal.Open(dolphin_vel_file).ReadAsArray()
            # Convert 0s to NaN
            vel_arr = np.where(vel_arr == 0, np.nan, vel_arr)
            # apply reference point
            if ref_meta is not None:
                ref_y = int(ref_meta['REF_Y'])
                ref_x = int(ref_meta['REF_X'])
                vel_arr -= np.nan_to_num(vel_arr[ref_y, ref_x])
            # write to file
            f["velocity"][:] = vel_arr

        print("finished writing to HDF5 file: {}".format(vel_file))

        recommended_mask_file = os.path.join(inps.out_dir, 'recommended_mask.h5')
        # generate mask file from unw phase field
        if inps.zero_mask is True:
            msk_file = os.path.join(os.path.dirname(ts_name),
                'combined_msk.h5')
            iargs = [recommended_mask_file, '-o', msk_file, '--nonzero']
            generate_mask.main(iargs)

            # mask TS file, since reference_point adds offset back in masked field
            iargs = [ts_name, '--mask', msk_file]
            mask.main(iargs)

            # mask velocity file
            iargs = [vel_file, '--mask', msk_file]
            mask.main(iargs)

    print(f"Total processing time: {(time.time() - start_time)/60:.2f} minutes")
    print("Done.")
    return


if __name__ == "__main__":
    main(sys.argv[1:])
