#!/usr/bin/env python
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Talib Oliver Cabrerra, Scott Staniewicz          #
############################################################

# Standard library imports
import argparse
import datetime
import glob
import itertools
import os
import sys
from pathlib import Path
from typing import Sequence

# Third-party imports
import affine
import h5py
import numpy as np
import pyproj
import rasterio
from osgeo import gdal
from rasterio import CRS
from rasterio.warp import reproject, Resampling
from tqdm import tqdm
import networkx as nx

from datetime import datetime as dt, timedelta
import rioxarray
import xarray
from osgeo import osr
import rasterio
from rasterio.warp import reproject, Resampling
from affine import Affine
from shapely import wkt
from shapely.geometry import Point
import random
import requests
import asf_search as asf
import matplotlib.pyplot as plt

import RAiDER
from RAiDER.models.hrrr import HRRR
from RAiDER.processWM import prepareWeatherModel
from RAiDER.delay import tropo_delay
from RAiDER.llreader import BoundingBox
from RAiDER.losreader import Raytracing
from eof.download import download_eofs

import warnings
warnings.filterwarnings("ignore")

# Add the src directory to sys.path
sys.path.append(str(Path(__file__).parent / 'src'))

# Local application/library-specific imports
from dem_stitcher.stitcher import stitch_dem
from mintpy.cli import (
    generate_mask,
    mask,
    temporal_average,
    timeseries2velocity
)
from mintpy.reference_point import reference_point_attribute
from mintpy.utils import arg_utils, ptime, readfile, writefile
from mintpy.utils import utils as ut
from mintpy.utils.utils0 import (
    azimuth2heading_angle,
    calc_azimuth_from_east_north_obs
)
from opera_utils import get_dates
from pst_dolphin_utils import (
    BackgroundRasterWriter,
    create_external_files,
    datetime_to_float,
    estimate_velocity,
    full_suffix,
    get_raster_bounds,
    get_raster_crs,
    get_raster_gt,
    get_raster_xysize,
    HDF5StackReader,
    load_gdal,
    process_blocks,
    warp_to_match
)
from tile_mate import get_raster_from_tiles
from tile_mate.stitcher import DATASET_SHORTNAMES

OPERA_DATASET_ROOT = './'

####################################################################################
EXAMPLE = """example:

  prep_mintpy.py
      -m pst_output/static_CSLCs/
      -c "pst_output/dolphin_output/stitched_interferograms/*.zeroed.cor.tif"
      -u "pst_output/dolphin_output/stitched_interferograms/*.unw.zeroed.tif"
      --geom-dir pst_output/dolphin_output/stitched_interferograms/geometry
      --ref-lalo '19.2485991551617 -155.32285148610057'
      -o mintpy_output

"""  # noqa: E501

# """
# Scott TODO:
# - UTM_ZONE, EPSG from the stitched IFG (it won't work to get a single GSLC burst)
# - pixel size is wrong since we're taking range/azimuth size, instead of geocoded size
# - HEIGHT: do we wanna try to get that from the saved orbit info?


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
        "-c",
        "--cor-file-glob",
        type=str,
        default="./interferograms/stitched/*.cor",
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
        "--water-mask-file",
        dest="water_mask_file",
        type=str,
        default=None,
        help="Specify either path to valid water mask, or download "
             f"using one of the following data sources: {DATASET_SHORTNAMES}",
    )
    parser.add_argument(
        "--dem-file",
        dest="dem_file",
        type=str,
        default=None,
        help="Specify either path to valid DEM, or download "
             "using one of the following data sources: srtm_v3, "
             "nasadem, glo_30, glo_90, glo_90_missing",
    )
    parser.add_argument(
        "--ref-lalo",
        dest="ref_lalo",
        type=str,
        default=None,
        help="Specify 'latitute longitude' of desired reference point. "
             "By default the pixel with the highest spatial coherence "
             "is selected",
    )
    parser.add_argument(
        "--min-coherence",
        dest="min_coherence",
        type=str,
        default='0.4',
        help="Specify minimum coherence of reference pixel "
             "for max-coherence method.",
    )
    parser.add_argument(
        "--zero-mask",
        dest="zero_mask",
        action="store_true",
        help="Mask all pixels with zero value in unw phase",
    )
    parser.add_argument(
        "--corr-lyrs",
        dest="corr_lyrs",
        action="store_true",
        help="Extract correction layers",
    )
    parser.add_argument(
        "--shortwvl-lyrs",
        dest="shortwvl_lyrs",
        action="store_true",
        help="Extract short wavelength layers",
    )

     # Add new arguments for tropospheric correction
    parser.add_argument(
        "--tropo-correction",
        dest="tropo_correction",
        action="store_true",
        help="Apply tropospheric correction using HRRR weather model"
    )
    parser.add_argument(
        "--work-dir",
        dest="work_dir",
        type=str,
        default="./work",
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

    t0 = datetime.datetime.strptime(
        meta_compass[f"{burst_ds}/sensing_start"][()].decode("utf-8"),
        "%Y-%m-%d %H:%M:%S.%f",
    )
    t1 = datetime.datetime.strptime(
        meta_compass[f"{burst_ds}/sensing_stop"][()].decode("utf-8"),
        "%Y-%m-%d %H:%M:%S.%f",
    )
    t_mid = t0 + (t1 - t0) / 2.0
    meta["CENTER_LINE_UTC"] = (
        t_mid - datetime.datetime(t_mid.year, t_mid.month, t_mid.day)
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


def write_coordinate_system(
    filename, dset_name, xy_dim_names=("x", "y"), grid_mapping_dset="spatial_ref"
):
    """Write the coordinate system CF metadata to an existing HDF5 file."""
    x_dim_name, y_dim_name = xy_dim_names
    atr = readfile.read_attribute(filename)
    epsg = int(atr.get("EPSG", 4326))

    with h5py.File(filename, "a") as hf:
        crs = pyproj.CRS.from_user_input(epsg)
        dset = hf[dset_name]

        # Setup the dataset holding the SRS information
        srs_dset = hf.require_dataset(grid_mapping_dset, shape=(), dtype=int)
        srs_dset.attrs.update(crs.to_cf())
        dset.attrs["grid_mapping"] = grid_mapping_dset

        if "date" in hf:
            date_arr = [
                datetime.datetime.strptime(ds, "%Y%m%d")
                for ds in hf["date"][()].astype(str)
            ]
            days_since = [(d - date_arr[0]).days for d in date_arr]
            dt_dim = hf.create_dataset("time", data=days_since)
            dt_dim.make_scale()
            cf_attrs = dict(
                units=f"days since {str(date_arr[0])}", calendar="proleptic_gregorian"
            )
            dt_dim.attrs.update(cf_attrs)
            dset.dims[0].attach_scale(dt_dim)
            dset.dims[0].label = "time"
        else:
            dt_dim = date_arr = None
            # If we want to do something other than time as a 3rd dimension...
            #  We'll need to figure out what other valid dims there are
            # otherwise, we can just do `phony_dims="sort"` in xarray

        # add metadata to x,y coordinates
        is_projected = crs.is_projected
        is_geographic = crs.is_geographic
        x_arr, y_arr = _get_xy_arrays(atr)
        x_dim_dset = hf.create_dataset(x_dim_name, data=x_arr)
        x_dim_dset.make_scale(x_dim_name)
        y_dim_dset = hf.create_dataset(y_dim_name, data=y_arr)
        y_dim_dset.make_scale(y_dim_name)

        x_coord_attrs = {}
        x_coord_attrs["axis"] = "X"
        y_coord_attrs = {}
        y_coord_attrs["axis"] = "Y"
        if is_projected:
            units = "meter"
            # X metadata
            x_coord_attrs["long_name"] = "x coordinate of projection"
            x_coord_attrs["standard_name"] = "projection_x_coordinate"
            x_coord_attrs["units"] = units
            # Y metadata
            y_coord_attrs["long_name"] = "y coordinate of projection"
            y_coord_attrs["standard_name"] = "projection_y_coordinate"
            y_coord_attrs["units"] = units
        elif is_geographic:
            # X metadata
            x_coord_attrs["long_name"] = "longitude"
            x_coord_attrs["standard_name"] = "longitude"
            x_coord_attrs["units"] = "degrees_east"
            # Y metadata
            y_coord_attrs["long_name"] = "latitude"
            y_coord_attrs["standard_name"] = "latitude"
            y_coord_attrs["units"] = "degrees_north"
        y_dim_dset.attrs.update(y_coord_attrs)
        x_dim_dset.attrs.update(x_coord_attrs)

        ndim = dset.ndim
        dset.dims[ndim - 1].attach_scale(x_dim_dset)
        dset.dims[ndim - 2].attach_scale(y_dim_dset)
        dset.dims[ndim - 1].label = x_dim_name
        dset.dims[ndim - 2].label = y_dim_name


def _get_date_pairs(filenames):
    str_list = [Path(f).stem for f in filenames]
    basenames_noext = [str(f).replace(full_suffix(f), "") for f in str_list]

    # access dates from beta products differently than from golden outputs
    date_pairs = []
    for i in basenames_noext:
        num_parts = i.split('_')
        if len(num_parts) == 9:
            date_pair = f'{num_parts[6][:8]}_{num_parts[7][:8]}'
            date_pairs.append(date_pair)
        if len(num_parts) == 2:
            date_pairs.append(i)

    return date_pairs


def get_azimuth_ang(dsDict):
    """Compute the azimuth angle from east/north coefficients"""
    east = dsDict["los_east"]
    north = dsDict["los_north"]
    azimuth_angle = calc_azimuth_from_east_north_obs(east, north)
    return azimuth_angle, east, north


def save_stack(
    fname,
    ds_name_dict,
    meta,
    file_list,
    water_mask,
    date12_list,
    phase2range=1,
    ref_y=None,
    ref_x=None,
    unw_file=False,
    mask_dict={},
    ):
    """Prepare h5 file for input stack of layers"""
    # initiate file
    writefile.layout_hdf5(fname, ds_name_dict, metadata=meta)

    # writing data to HDF5 file
    reflyr_name = file_list[0].split(':')[-1]
    print("writing data to HDF5 file {} with a mode ...".format( \
          fname))
    num_file = len(file_list)
    with h5py.File(fname, "a") as f:
        prog_bar = ptime.progressBar(maxValue=num_file)
        for i,file_inc in enumerate(file_list):
            # read data using gdal
            data = load_gdal(file_inc, masked=True)

            # apply reference point, if not None
            if ref_y is not None and ref_x is not None:
                data -= np.nan_to_num(data[ref_y, ref_x])

            # mask by specified dict of thresholds
            for dict_key in mask_dict.keys():
                mask_lyr = file_inc.replace(reflyr_name, dict_key)
                mask_thres = mask_dict[dict_key]
                mask_data = load_gdal(mask_lyr)
                data[mask_data == mask_thres] = np.nan

            # also apply water mask
            data = data * water_mask

            # if unw file, convert nans to 0
            # necessary to avoid errors with MintPy velocity fitting
            if unw_file is True:
                data = np.nan_to_num(data)

            # apply conversion factor and add to cube
            f["timeseries"][i + 1] = data * phase2range
            prog_bar.update(i + 1, suffix=date12_list[i])

        prog_bar.close()

        print("set value at the first acquisition to ZERO.")
        f["timeseries"][0] = 0.0

    print("finished writing to HDF5 file: {}".format(fname))

def prepare_timeseries(
    outfile,
    unw_files,
    track_version,
    metadata,
    water_mask_file=None,
    ref_lalo=None,
    corr_lyrs=False,
    shortwvl_lyrs=False,
    apply_tropo_correction=False,
    median_height=50,
    work_dir=None
):
    """
    Prepare the timeseries file accounting for different reference dates
    in input files.
    """
    print("-" * 50)
    print("preparing timeseries file: {}".format(outfile))

    # copy metadata to meta
    meta = {key: value for key, value in metadata.items()}

    # pass lyr to disp conversion factor, which depends on the product version
    sp_coh_lyr_name = 'interferometric_correlation'
    if track_version == 0.3:
        disp_lyr_name = 'unwrapped_phase'
        phase2range = -1 * float(meta["WAVELENGTH"]) / (4.0 * np.pi)
    if track_version >= 0.4:
        disp_lyr_name = 'displacement'
        phase2range = 1
    if track_version == 0.7:
       sp_coh_lyr_name = 'estimated_spatial_coherence'
    if track_version == 0.8:
       sp_coh_lyr_name = 'estimated_phase_quality'

    if apply_tropo_correction and work_dir:
        os.makedirs(work_dir, exist_ok=True)
        os.makedirs(f'{work_dir}/orbits', exist_ok=True)
        os.makedirs(f'weather_files', exist_ok=True)

    # grab date list from the filename
    date12_list = _get_date_pairs(unw_files)
    num_file = len(unw_files)
    print("number of unwrapped interferograms: {}".format(num_file))

    # Create a directed graph to represent date connections
    G = nx.DiGraph()
    
    # Add edges (measurements) to the graph
    date_pairs = [dl.split("_") for dl in date12_list]
    for (ref_date, sec_date), file in zip(date_pairs, unw_files):
        G.add_edge(ref_date, sec_date, file=file)

    # Get all unique dates
    date_list = sorted(set(itertools.chain.from_iterable(date_pairs)))
    num_date = len(date_list)
    print("number of acquisitions: {}\n{}".format(num_date, date_list))

    # size info
    cols, rows = get_raster_xysize(unw_files[0])

    # baseline info
    pbase = np.zeros(num_date, dtype=np.float32)

    # define dataset structure
    dates = np.array(date_list, dtype=np.string_)
    ds_name_dict = {
        "date": [dates.dtype, (num_date,), dates],
        "bperp": [np.float32, (num_date,), pbase],
        "timeseries": [np.float32, (num_date, rows, cols), None],
    }

    # read water mask
    if water_mask_file is not None:
         water_mask = readfile.read(water_mask_file,
             datasetName='waterMask')[0]
    else:
        water_mask = np.ones((rows, cols), dtype=np.float32)

    # handle reference point
    coord = ut.coordinate(meta)
    if ref_lalo is not None:
        ref_lat, ref_lon = ref_lalo.split()
        ref_y, ref_x = coord.geo2radar(np.array(float(ref_lat)),
                                     np.array(float(ref_lon)))[0:2]
    else:
        coh_file = os.path.join(os.path.dirname(outfile), 'avgSpatialCoh.h5')
        coh = readfile.read(coh_file)[0]
        if water_mask.dtype.name == 'bool':
            coh[water_mask == False] = 0.0
        else:
            coh[water_mask == 0] = 0.0
        ref_y, ref_x = np.unravel_index(np.argmax(coh), coh.shape)
        del coh

    # update metadata
    ref_meta = reference_point_attribute(meta, y=ref_y, x=ref_x)
    meta.update(ref_meta)
    meta["FILE_TYPE"] = "timeseries"
    meta["UNIT"] = "m"

    def calculate_cumulative_displacement(
        date, water_mask, mask_dict, reflyr_name
    ):
        """
        Calculate cumulative displacement up to given date using shortest path
        """
        if date == date_list[0]:
            return np.zeros((rows, cols), dtype=np.float32)
        
        # Find shortest path from first date to current date
        try:
            # Finding all paths (files) in the shortest path from
            # the first date to given date
            path = nx.shortest_path(G, source=date_list[0], target=date)
        except nx.NetworkXNoPath:
            print(f"Warning: No path found to date {date}")
            return None

        # Calculate cumulative displacement along the path
        cumulative = np.zeros((rows, cols), dtype=np.float32)
        print(f'\ndate for calculating cumulative displacement: {date}')
        print(f'shortest path from {date_list[0]} to {date}')
        for i in range(len(path)-1):
            ref_date, sec_date = path[i], path[i+1]
            file = G[ref_date][sec_date]['file']	# File in the shortest path
            print(f'{i} {ref_date} to {sec_date}: {file}')	

            # Read displacement data
            data = load_gdal(file, masked=True)
            data *= phase2range

            # Apply tropospheric correction if requested
            if apply_tropo_correction and work_dir:
                print(f"\nApplying tropospheric correction to pair {date12_list[i]}")
                
                # Extract original NetCDF filename from the GDAL path string
                # Format is 'NETCDF:"path/to/file.nc":displacement'
                nc_file = file.split('"')[1]
                
                try:
                    # Read parameters needed for tropo correction
                    with h5py.File(nc_file, 'r') as nc:
                        track_number = nc['identification']['track_number'][()]
                        bounding_polygon = nc['identification']['bounding_polygon'][()].decode()
                        ref_datetime = dt.strptime(nc['identification']['reference_datetime'][()].decode(),
                                                 '%Y-%m-%d %H:%M:%S.%f')
                        sec_datetime = dt.strptime(nc['identification']['secondary_datetime'][()].decode(),
                                                 '%Y-%m-%d %H:%M:%S.%f')
                        spatial_ref_attrs = nc['spatial_ref'].attrs
                        crs_wkt = spatial_ref_attrs['crs_wkt']   
                        epsg_code = crs_wkt.split('ID["EPSG",')[-1].split(']')[0]
                        epsg_str = f'EPSG:{epsg_code}'
                        GeoTransform = spatial_ref_attrs['GeoTransform']
                        frame_id = nc['identification']['frame_id'][()]
                        frame_id = 'F' + str(frame_id).zfill(5)
                        mission_id = nc['identification']['mission_id'][()]
                        ref_date = ref_datetime.strftime('%Y%m%d')  # YYYYMMDD
                        sec_date = sec_datetime.strftime('%Y%m%d')
                        unwrapper_mask = nc['unwrapper_mask'][:]

                    # Setup parameters for tropospheric correction
                    params = {
                        'track_number': track_number,
                        'bounding_polygon': bounding_polygon,
                        'ref_datetime': ref_datetime,
                        'sec_datetime': sec_datetime,
                        'GeoTransform': GeoTransform,
                        'epsg' : epsg_str,
                        'median_height' : median_height,
                        'mission_id' : mission_id,
                        'height' : data.shape[0],
                        'width' : data.shape[1],
                    }

                    # Calculate and apply tropospheric correction
                    calculated_tropo_delay = calculate_tropospheric_delay(params, work_dir)     # unit: meter
                    
                    ### code to plot the tropospheric correction
                    _data = data.copy()
                    _data[unwrapper_mask==0.] = np.nan
                    fig, ax = plt.subplots(1, 3, figsize=(15,10))
                    im0 = ax[0].imshow(_data, cmap='RdBu')
                    ax[0].set_title(f'Displacement before tropo correction \n{frame_id} {ref_date}-{sec_date}')
                    ax[0].axis('off')
                    plt.colorbar(im0, ax=ax[0], label='LOS (m)', shrink=0.2)
                    im1 = ax[1].imshow(calculated_tropo_delay, cmap='RdBu')
                    ax[1].set_title(f'Tropospheric delay \n{frame_id} {ref_date}-{sec_date}')
                    ax[1].axis('off')
                    plt.colorbar(im1, ax=ax[1], label='LOS (m)',shrink=0.2)
                    im2 = ax[2].imshow(_data - calculated_tropo_delay, cmap='RdBu')
                    ax[2].set_title(f'Displacement after tropo correction \n{frame_id} {ref_date}-{sec_date}')
                    ax[2].axis('off')
                    plt.colorbar(im2, ax=ax[2], label='LOS (m)', shrink=0.2)
                    plt.tight_layout()
                    fig.savefig(f'{work_dir}/tropo_corrected_displacement_{frame_id}_Raytracing_{ref_date}_{sec_date}.png', dpi=300, bbox_inches='tight')
                    plt.close('all')
                    del _data
                    del unwrapper_mask
                    ###

                    data -= calculated_tropo_delay      # unit: meter 
 
                except Exception as e:
                    print(f"Warning: Tropospheric correction failed for {ref_date}_{sec_date}: {str(e)}")
                    print("Continuing with uncorrected data...")

            # Apply reference point correction
            if ref_y is not None and ref_x is not None:
                data -= np.nan_to_num(data[ref_y, ref_x])
          
            # mask by specified dict of thresholds 
            for dict_key in mask_dict.keys():
                mask_lyr = file.replace(reflyr_name, dict_key)
                mask_thres = mask_dict[dict_key]
                mask_data = load_gdal(mask_lyr)
                data[mask_data < mask_thres] = np.nan
 
            # Apply water mask
            data *= water_mask

            # Add to cumulative displacement
            cumulative += data 

        # convert nans to 0
        # necessary to avoid errors with MintPy velocity fitting
        return np.nan_to_num(cumulative)

    # Write timeseries data
    print(f"Writing data to HDF5 file {outfile}")
    writefile.layout_hdf5(outfile, ds_name_dict, metadata=meta)

    # set dictionary that will be used to mask TS by specified thresholds
    mask_dict = {}
    mask_dict['connected_component_labels'] = 1
    mask_dict['temporal_coherence'] = 0.6
    mask_dict[sp_coh_lyr_name] = 0.5
    if track_version == 0.8:
        mask_dict['water_mask'] = 1

    reflyr_name = unw_files[0].split(':')[-1]
    
    with h5py.File(outfile, "a") as f:
        prog_bar = ptime.progressBar(maxValue=num_date)
        
        # First date is always zero
        f["timeseries"][0] = np.zeros((rows, cols), dtype=np.float32)
        
        # Calculate cumulative displacement for each date
        for i, date in enumerate(date_list[1:], start=1):
            displacement = calculate_cumulative_displacement(
                date, water_mask, mask_dict, reflyr_name)
            if displacement is not None:
                # writing timeseries to the date
                f["timeseries"][i] = displacement
            prog_bar.update(i, suffix=date)
        
        prog_bar.close()

    print("finished writing to HDF5 file: {}".format(outfile))
    
    # Handle additional layers (correction layers, mask layers, etc.)
    all_outputs = [outfile]
    mask_layers = ['connected_component_labels', 'temporal_coherence',
                  sp_coh_lyr_name]

    correction_layers = []
    if corr_lyrs is True:
        correction_layers = ['ionospheric_delay',
                           'solid_earth_tide']
        if track_version < 0.8:
            correction_layers.append('tropospheric_delay')

    # Process additional layers similarly to main displacement
    for lyr in mask_layers + correction_layers:
        lyr_fname = os.path.join(os.path.dirname(outfile), f'{lyr}.h5')
        if lyr in correction_layers:
            lyr_paths = [i.replace(disp_lyr_name, f'/corrections/{lyr}')
                        for i in unw_files]
            save_stack(lyr_fname, ds_name_dict, meta, lyr_paths,
                      water_mask, date12_list, phase2range, ref_y, ref_x)
        else:
            lyr_paths = [i.replace(disp_lyr_name, lyr) for i in unw_files]
            save_stack(lyr_fname, ds_name_dict, meta, lyr_paths,
                      water_mask, date12_list, 1)
        all_outputs.append(lyr_fname)

    # Handle short wavelength layers
    if shortwvl_lyrs is True and track_version >= 0.4:
        shortwvl_lyr = 'short_wavelength_displacement'
        shortwvl_fname = os.path.join(os.path.dirname(outfile),
            f'{shortwvl_lyr}.h5')
        shortwvl_files = [i.replace(disp_lyr_name, shortwvl_lyr)
                         for i in unw_files]
        save_stack(shortwvl_fname, ds_name_dict, meta,
            shortwvl_files, water_mask, date12_list, phase2range,
            ref_y, ref_x)
        all_outputs.append(shortwvl_fname)

    return all_outputs, ref_meta

def mintpy_prepare_geometry(outfile, geom_dir, metadata,
                            water_mask_file=None):
    """Prepare the geometry file."""
    print('-' * 50)
    print(f'preparing geometry file: {outfile}')

    geom_path = Path(geom_dir)
    # copy metadata to meta
    meta = {key: value for key, value in metadata.items()}
    meta["FILE_TYPE"] = "geometry"

    file_to_path = {
        'los_east': geom_path / 'los_east.tif',
        'los_north': geom_path / 'los_north.tif',
        'height': geom_path / 'height.tif',
        'shadowMask': geom_path / 'layover_shadow_mask.tif',
    }

    if water_mask_file:
        file_to_path['waterMask'] = water_mask_file

    dsDict = {}
    for dsName, fname in file_to_path.items():
        try:
            data = readfile.read(fname, datasetName=dsName)[0]
            # TODO: add general functionality to handle nodata into Mintpy
            if dsName not in ['shadowMask', 'waterMask']:
                data[data == 0] = np.nan
            dsDict[dsName] = data

            # write data to HDF5 file
        except FileNotFoundError as e:  # https://github.com/insarlab/MintPy/issues/1081
            print(f'Skipping {fname}: {e}')

    # Compute the azimuth and incidence angles from east/north coefficients
    azimuth_angle, east, north = get_azimuth_ang(dsDict)
    dsDict['azimuthAngle'] = azimuth_angle

    up = np.sqrt(1 - east**2 - north**2)
    incidence_angle = np.rad2deg(np.arccos(up))
    dsDict['incidenceAngle'] = incidence_angle

    writefile.write(dsDict, outfile, metadata=meta)
    return outfile

def chunked_nanmean(file_list, mask, chunk_size=50):
    ''' averaging over chunks for handling memory efficiently '''
    dat = load_gdal(file_list[0], masked=True) * mask
    row, col = dat.shape  
   
    total_sum = np.zeros((row,col))
    total_count = np.zeros((row,col))  

    # process files in chunks
    for i in tqdm(range(0, len(file_list), chunk_size)):
        chunk_files = file_list[i:i+chunk_size]
        chunk_data = []

        # read data from each file in the chunk
        for file in chunk_files:
            chunk_data.append(load_gdal(file, masked=True) * mask)

        # stack the chunk data
        chunk = np.stack(chunk_data)
      
        # calculate sum and count for the chunk 
        chunk_sum = np.nansum(chunk, axis=0)
        chunk_count = np.sum(~np.isnan(chunk), axis=0)  

        # update total sum and count
        total_sum += chunk_sum
        total_count += chunk_count

    # calculate the mean
    result = total_sum / total_count

    return result
         
def prepare_average_stack(outfile, infiles, water_mask, file_type, metadata):
    """Average and export specified layers"""
    print("-" * 50)
    print("preparing average of stack: {}".format(outfile))

    # copy metadata to meta
    meta = {key: value for key, value in metadata.items()}
    meta["FILE_TYPE"] = file_type # "mask" "temporalCoherence"
    meta["UNIT"] = "1"

    # calculating nanmean in chunks
    data = chunked_nanmean(infiles, water_mask)

    # write to HDF5 file
    writefile.write(data, outfile, metadata=meta)
    return outfile

def calculate_tropospheric_delay(params, work_dir):
    """Calculate tropospheric delay for a single interferogram"""
    
    # Initialize weather model
    weather_model = HRRR()
    
    # Extract parameters
    track_number = params['track_number']
    bounding_polygon = params['bounding_polygon']
    ref_datetime = params['ref_datetime']
    sec_datetime = params['sec_datetime']
    GeoTransform = params['GeoTransform']
    epsg_str = params['epsg'] 
    mission_id = params['mission_id']
    height = params['height']
    width = params['width']
    evaluation_height = int(params['median_height'])    # Get evaluation height

    # Create polygon object
    polygon = wkt.loads(bounding_polygon)
    datetime_list = [ref_datetime, sec_datetime]

    transform = [float(x) for x in GeoTransform.split()]
    x_res = transform[1] 
    y_res = transform[5] 
    x_min = transform[0]  
    y_max = transform[3]  
 
    corners = get_bounds_from_geotransform(GeoTransform, width, height)     # bounds of input DISP-S1 based on GeoTransform
    transformed_corners = transform_coords(corners, epsg_str, 'EPSG:4326')  # convert to EPSG:4326

    # Calculate bounds
    lons, lats = zip(*transformed_corners)
    min_lat, max_lat = min(lats), max(lats)
    min_lon, max_lon = min(lons), max(lons)

    ll_bounds = [min_lat, max_lat, min_lon, max_lon]   # bounds of the DISP-S1  
    print('Bounds of DISP-S1', ll_bounds)

    # Process each date
    delays = []
    for _datetime in datetime_list:

        _date = _datetime.strftime('%Y%m%d')
        
        # Find sensor
        sensor = find_sentinel1_sensor(bounding_polygon, int(track_number), _date)
        if not sensor:
            sensor = mission_id
            print('WARNING: sensor name was not correctly found by asf_search\n')
            # raise ValueError(f"Could not determine sensor for date {_date}")
            
        # Setup area of interest
        aoi = BoundingBox(ll_bounds)
        aoi.add_buffer(weather_model.getLLRes())
        aoi.set_output_xygrid(4326)
        
        # Download orbit file
        orbit_file = str(download_eofs([_datetime], [sensor], 
                                     save_dir=os.path.join(work_dir, 'orbits'),
                                     orbit_type='precise')[0])
        
        print('Downloaded orbit file: ', orbit_file)
        
        # Setup LOS calculator
        los = Raytracing(orbit_file, time=_datetime)
        
        if los.ray_trace():
            wm_bounds = aoi.calc_buffer_ray(
                los.getSensorDirection(),
                lookDir=los.getLookDirection(),
                incAngle=30
            )
        else:
            wm_bounds = aoi.bounds()
            
        # Process weather model
        weather_nc = prepareWeatherModel(
            weather_model, 
            _datetime,
            ll_bounds=wm_bounds,
            makePlots=False
        )
        
        weather_model.set_latlon_bounds(
            wm_bounds,
            output_spacing=aoi.get_output_spacing()
        )
        
        # Calculate delays
        ds, _ = tropo_delay(
            _datetime,
            weather_model.out_file('weather_files'),
            aoi,
            los,
            height_levels=[0, 100, 500, 1000]
        )

        delays.append((ds['wet'] + ds['hydro']).interp(z=evaluation_height))
    
    # Calculate differential delay
    diff_delays_ds = (delays[0] - delays[1])    # If RAiDER output is two-way delay 
    # diff_delays_ds = (delays[0] - delays[1]) * 2   # If RAiDER output is one-way delay
    
    # Rename the coordinates to follow CF conventions
    diff_delays_ds = diff_delays_ds.rename({'x': 'longitude', 'y': 'latitude'})

    # Add proper attributes
    diff_delays_ds['longitude'].attrs.update({
        'units': 'degrees_east',
        'standard_name': 'longitude',
        'long_name': 'longitude'
    })

    diff_delays_ds['latitude'].attrs.update({
        'units': 'degrees_north',
        'standard_name': 'latitude',
        'long_name': 'latitude'
    })

    # Add CRS information
    diff_delays_ds.attrs['crs'] = 'EPSG:4326'
    
    # Export and reproject
    temp_tif = os.path.join(work_dir, 'temp_delay.tif')
    diff_delays_ds.rio.to_raster(temp_tif)
    
    with rasterio.open(temp_tif) as src:
        transform = Affine(x_res, 0.0, x_min,
                      0.0, y_res, y_max)
        
        meta = src.meta.copy()
        meta.update({
            'crs': epsg_str,
            'transform': transform,
            'width': width,
            'height': height
        })
        
        reprojected_tif = os.path.join(work_dir, 'reprojected_delay.tif')
        with rasterio.open(reprojected_tif, 'w', **meta) as dst:
            reproject(
                source=rasterio.band(src, 1),
                destination=rasterio.band(dst, 1),
                src_transform=src.transform,
                src_crs=src.crs,
                dst_transform=transform,
                dst_crs=meta['crs'],
                resampling=Resampling.lanczos
            )
    
    # Read reprojected delay
    with rasterio.open(reprojected_tif) as src:
        tropo_delay_arr = src.read(1)
    
    # Clean up temporary files
    os.remove(temp_tif)
    os.remove(reprojected_tif)
    
    return tropo_delay_arr

def find_sentinel1_sensor(wkt_polygon, track_number, start_date):
    """
    Search for Sentinel-1 scenes and return the sensor (S1A or S1B) of the first found scene.
    
    Args:
        wkt_polygon (str): WKT representation of the area of interest
        track_number (int): track number of Sentinel-1
        start_date (str): Start date in 'YYYYMMDD' format
    
    Returns:
        str: Sensor name ('S1A' or 'S1B') or None if no scenes found
    """
    try:
        # Convert dates to datetime objects
        start = dt.strptime(start_date, '%Y%m%d')
        end = start + timedelta(days=1)

        results = asf.search(
            platform=[asf.PLATFORM.SENTINEL1],
            relativeOrbit=[track_number],
            processingLevel=[asf.PRODUCT_TYPE.SLC],
            start=start,
            end=end,
            maxResults=5,
            intersectsWith=wkt_polygon
        )

        if len(results) > 0:
            # Get the first scene's properties
            first_scene = results[0]
            
            # Extract platform name (Sentinel-1A or Sentinel-1B)
            sensor = first_scene.properties['platform']

            sensor_mapping = {
                "Sentinel-1A": "S1A",
                "Sentinel-1B": "S1B"
            }

            return sensor_mapping[sensor]
        else:
            print("No Sentinel-1 scenes found for the specified criteria")
            return None
            
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        return None

def get_bounds_from_geotransform(geotransform, width, height):
    """
    Get the corner coordinates from a GeoTransform
    """
    # GeoTransform format: (top_left_x, pixel_width, x_rotation, top_left_y, y_rotation, pixel_height)
    gt = [float(x) for x in geotransform.strip("'").split()]
    
    # Calculate corner coordinates in the original projection
    corners = [
        (gt[0], gt[3]),  # Top-left
        (gt[0] + width * gt[1], gt[3]),  # Top-right
        (gt[0], gt[3] + height * gt[5]),  # Bottom-left
        (gt[0] + width * gt[1], gt[3] + height * gt[5])  # Bottom-right
    ]
    return corners

def transform_coords(coords, source_epsg, target_epsg):
    """
    Transform coordinates from source EPSG to target EPSG
    """
    # Create coordinate transformation
    source = osr.SpatialReference()
    source.ImportFromEPSG(int(source_epsg.split(':')[1]))
    
    target = osr.SpatialReference()
    target.ImportFromEPSG(int(target_epsg.split(':')[1]))
    
    transform = osr.CoordinateTransformation(source, target)
    
    # Transform all coordinates
    transformed_coords = []
    for x, y in coords:
        point = transform.TransformPoint(x, y)
        transformed_coords.append((point[1], point[0])) 
    return transformed_coords

def prepare_stack(
    outfile,
    unw_files,
    disp_lyr_name,
    track_version,
    metadata,
    water_mask_file=None,
    ref_lalo=None,
):
    """Prepare the input unw stack with tropospheric correction."""
    print("-" * 50)
    print("preparing ifgramStack file: {}".format(outfile))
    
    # copy metadata to meta
    meta = {key: value for key, value in metadata.items()}

    # get list of *.unw file
    num_pair = len(unw_files)

    conv_factor = 1
    # >=v0.4 increment in units of m, must be converted to phs
    sp_coh_lyr_name = 'interferometric_correlation'
    if track_version >= 0.4:
        conv_factor = -1 * (4.0 * np.pi) / float(metadata["WAVELENGTH"])
    if track_version == 0.7:
       sp_coh_lyr_name = 'estimated_spatial_coherence'
    if track_version == 0.8:
       sp_coh_lyr_name = 'estimated_phase_quality'

    print(f"number of unwrapped interferograms: {num_pair}")

    # get list of interferometric correlation layers
    cor_files = [i.replace(disp_lyr_name, sp_coh_lyr_name) for i in unw_files]
    print(f"number of correlation files: {len(cor_files)}")

    # get list of conncomp layers
    cc_files = [i.replace(disp_lyr_name, 'connected_component_labels') for i in unw_files]
    print(f"number of connected components files: {len(cc_files)}")

    if len(cc_files) != len(unw_files) or len(cor_files) != len(unw_files):
        print("the number of *.unw and *.unw.conncomp or *.cor files are NOT consistent")
        if len(unw_files) > len(cor_files):
            print("skip creating ifgramStack.h5 file.")
            return

        print("Keeping only cor files which match a unw file")
        unw_dates_set = set([tuple(get_dates(f)) for f in unw_files])
        cor_files = [f for f in cor_files if tuple(get_dates(f)) in unw_dates_set]

    # get date info: date12_list
    date12_list = _get_date_pairs(unw_files)

    # TODO: compute the spatial baseline using COMPASS metadata
    pbase = np.zeros(num_pair, dtype=np.float32)

    # size info
    cols, rows = get_raster_xysize(unw_files[0])

    # define (and fill out some) dataset structure
    date12_arr = np.array([x.split("_") for x in date12_list], dtype=np.string_)
    drop_ifgram = np.ones(num_pair, dtype=np.bool_)
    ds_name_dict = {
        "date": [date12_arr.dtype, (num_pair, 2), date12_arr],
        "bperp": [np.float32, (num_pair,), pbase],
        "dropIfgram": [np.bool_, (num_pair,), drop_ifgram],
        "unwrapPhase": [np.float32, (num_pair, rows, cols), None],
        "coherence": [np.float32, (num_pair, rows, cols), None],
        "connectComponent": [np.float32, (num_pair, rows, cols), None],
    }

    # read water mask
    if water_mask_file is not None:
         water_mask = readfile.read(water_mask_file, datasetName='waterMask')[0]
    else:
        water_mask = np.ones((rows, cols), dtype=np.float32)

    # initiate HDF5 file
    meta["FILE_TYPE"] = "ifgramStack"
    writefile.layout_hdf5(outfile, ds_name_dict, metadata=meta)

    # writing data to HDF5 file
    print("writing data to HDF5 file {} with a mode ...".format(outfile))
    with h5py.File(outfile, "a") as f:
        prog_bar = ptime.progressBar(maxValue=num_pair)
        for i, (unw_file, cor_file, cc_file) in enumerate(zip(unw_files, cor_files, cc_files)):
            # read/write *.unw file
            unw_arr = load_gdal(unw_file, masked=True) * water_mask * conv_factor 
            # mask
            unw_arr[unw_arr == 0] = np.nan

            f["unwrapPhase"][i] = unw_arr     # unit: radian

            # read/write *.cor file
            f["coherence"][i] = load_gdal(cor_file, masked=True) * water_mask

            # read/write *.unw.conncomp file
            f["connectComponent"][i] = load_gdal(cc_file, masked=True) * water_mask

            prog_bar.update(i + 1, suffix=date12_list[i])
        prog_bar.close()

    print("finished writing to HDF5 file: {}".format(outfile))

    # generate average spatial coherence
    coh_dir = os.path.dirname(os.path.dirname(outfile))
    coh_file = os.path.join(coh_dir, 'avgSpatialCoh.h5')
    iargs = [outfile, '--dataset', 'coherence', '-o', coh_file, '--update']
    temporal_average.main(iargs)

    if water_mask_file is not None:
        # determine whether the defined reference point is valid
        if ref_lalo is not None:
            ref_lat, ref_lon = ref_lalo.split()
            # get value at coordinate
            atr = readfile.read_attribute(coh_file)
            coord = ut.coordinate(atr)
            (ref_y,
             ref_x) = coord.geo2radar(np.array(float(ref_lat)),
                                      np.array(float(ref_lon)))[0:2]
            val_at_refpoint = water_mask[ref_y, ref_x]
            # exit if reference point is in masked area
            if val_at_refpoint == False or val_at_refpoint ==  0:
                raise Exception(f'Specified input --ref-lalo {ref_lalo} '
                                'not in masked region. Inspect output file'
                                f'{water_mask_file} to inform selection '
                                'of new point.')

    # extract correction layers
    # Remove dict entries not relevant to correction layers
    for key in ['coherence', 'connectComponent']:
        ds_name_dict.pop(key, None)

    # loop through and create files for persistent scatterer
    ps_files = \
            [i.replace(disp_lyr_name, 'persistent_scatterer_mask') \
             for i in unw_files]
    ps_fname = os.path.join(os.path.dirname(outfile), 
        'persistent_scatterer_mask.h5')
    prepare_average_stack(ps_fname, ps_files, water_mask, 'mask', meta)

    # loop through and create files for temporal coherence
    tempcoh_files = \
            [i.replace(disp_lyr_name, 'temporal_coherence') \
             for i in unw_files]
    tempcoh_fname = os.path.join(os.path.dirname(outfile),
        'temporalCoherence.h5')
    prepare_average_stack(tempcoh_fname, tempcoh_files, water_mask,
                          'temporalCoherence', meta)

    # loop through and create files for temporal coherence
    conn_files = \
            [i.replace(disp_lyr_name, 'connected_component_labels') \
             for i in unw_files]
    conn_fname = os.path.join(os.path.dirname(outfile),
        'connectedComponent.h5')
    prepare_average_stack(conn_fname, conn_files, water_mask,
                          'connectedComponent', meta)

    # loop through and create files for temporal coherence
    spcoh_files = \
            [i.replace(disp_lyr_name, sp_coh_lyr_name) \
             for i in unw_files]
    spcoh_fname = os.path.join(os.path.dirname(outfile),
        'estimatedSpatialCoherence.h5')
    prepare_average_stack(spcoh_fname, spcoh_files, water_mask,
                          'estimatedSpatialCoherence', meta)

    return

def main(iargs=None):
    """Run the preparation functions."""
    inps = cmd_line_parse(iargs)

    unw_files = sorted(
        glob.glob(inps.unw_file_glob),
        key=lambda x: x.split('_')[7][:8]
    )

    # filter input by specified dates
    if inps.startDate is not None:
        startDate = int(inps.startDate)
        filtered_unw_files = []
        for i in unw_files:
            sec_date = int(os.path.basename(i).split('_')[7][:8])
            if sec_date >= startDate:
                filtered_unw_files.append(i)
        unw_files = filtered_unw_files

    if inps.endDate is not None:
        endDate = int(inps.endDate)
        filtered_unw_files = []
        for i in unw_files:
            sec_date = int(os.path.basename(i).split('_')[7][:8])
            if sec_date <= endDate:
                filtered_unw_files.append(i)
        unw_files = filtered_unw_files
    
    date12_list = _get_date_pairs(unw_files)
    print(f"Found {len(unw_files)} unwrapped files")

    # track product version
    track_version = []
    for i in unw_files:
        fname = os.path.basename(i)
        version_n = float(fname.split('_')[-2].split('v')[1])
        # round to nearest tenth
        version_n = round(version_n, 1)
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
    if track_version == 0.3:
        disp_lyr_name = 'unwrapped_phase'
    if track_version >= 0.4:
        disp_lyr_name = 'displacement'

    # append appropriate NETCDF prefixes
    unw_files = \
        [f'NETCDF:"{i}":{disp_lyr_name}' for i in unw_files]

    # get geolocation info
    static_dir = Path(inps.meta_file)
    geometry_dir = Path(inps.geom_dir)
    geometry_dir.mkdir(exist_ok=True)
    crs = get_raster_crs(unw_files[0])
    epsg = crs.to_epsg()
    out_bounds = get_raster_bounds(unw_files[0])

    # create water mask, if not specified
    Path(inps.out_dir).mkdir(exist_ok=True)
    if inps.water_mask_file is not None:
        if not Path(inps.water_mask_file).exists():
            inps.water_mask_file = create_external_files(inps.water_mask_file,
               unw_files[0], out_bounds, crs, inps.out_dir,
               maskfile=True)

    # create DEM, if not specified
    if inps.dem_file is not None:
        if not Path(inps.dem_file).exists():
            inps.dem_file = create_external_files(inps.dem_file,
               unw_files[0], out_bounds, crs, inps.out_dir,
               demfile=True)
    # check if mask file already generated through dolphin
    else:
        dem_file = geometry_dir.glob('*_DEM.tif')
        dem_file = [i for i in dem_file]
        if len(dem_file) > 0:
            inps.dem_file = dem_file[0]
            print(f"Found DEM file {inps.dem_file}")

    # check if geometry files exist
    geometry_files = ['los_east.tif', 'los_north.tif',
                      'layover_shadow_mask.tif']
    missing_files = [file for file in geometry_files
                     if not (geometry_dir / file).is_file()]
    if missing_files != []:
        missing_files_str = ', '.join(missing_files)
        raise FileNotFoundError(
            'The following geometry files are missing in the directory '
            f'{geometry_dir}: {missing_files_str}'
        )

    # check if height file exists
    height_file = geometry_dir / "height.tif"
    if not height_file.is_file():
        print(f"Generated {height_file} file")
        warp_to_match(
                    input_file=inps.dem_file,
                    match_file=unw_files[0],
                    output_file=height_file,
                    resample_alg="cubic",
                )
    
    print('DEM file: ', inps.dem_file)
    with rasterio.open(inps.dem_file) as src:
        height_arr = np.nan_to_num(src.read(1))
        valid_mask = (height_arr != 0)
        median_height = int(np.median(height_arr[valid_mask]))

    # check static layer naming convention
    allcaps_geometry = True
    static_files = sorted(Path(static_dir).glob("*STATIC_*.h5"))
    # capture alternate filename convention
    if static_files == []:
        allcaps_geometry = False
        static_files = sorted(Path(static_dir).glob("*static_*.h5"))

    # translate input options
    processor = "sweets"  # isce_utils.get_processor(inps.meta_file)
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

    # output directory
    for dname in [inps.out_dir, os.path.join(inps.out_dir, "inputs")]:
        os.makedirs(dname, exist_ok=True)

    # prepare ifgstack, which contains UNW phase + conn comp + coherence
    stack_file = os.path.join(inps.out_dir, "inputs/ifgramStack.h5")
    prepare_stack(
        outfile=stack_file,
        unw_files=unw_files,
        disp_lyr_name=disp_lyr_name,
        track_version=track_version,
        metadata=meta,
        water_mask_file=inps.water_mask_file,
        ref_lalo=inps.ref_lalo,
    )

    # prepare TS file
    og_ts_file = os.path.join(inps.out_dir, "timeseries.h5")
    geom_file = os.path.join(inps.out_dir, "geometryGeo.h5")
    all_outputs = [og_ts_file]
    ref_meta = None
    # time-series (if inputs are multi-reference and outputs are for single-reference)
    all_outputs, ref_meta = prepare_timeseries(
        outfile=og_ts_file,
        unw_files=unw_files,
        track_version=track_version,
        metadata=meta,
        water_mask_file=inps.water_mask_file,
        ref_lalo=inps.ref_lalo,
        corr_lyrs=inps.corr_lyrs,
        shortwvl_lyrs=inps.shortwvl_lyrs,
        apply_tropo_correction=inps.tropo_correction,
        median_height=median_height,
        work_dir=inps.work_dir if inps.tropo_correction else None
    )

    mintpy_prepare_geometry(geom_file, geom_dir=inps.geom_dir, metadata=meta,
                            water_mask_file=inps.water_mask_file)

    # generate velocity fit
    # first set variables
    dolphin_ref_tif = os.path.join(inps.out_dir, 'dolphin_reference.tif')
    dolphin_vel_file = os.path.join(inps.out_dir, 'velocity.tif')
    vel_file = os.path.join(inps.out_dir, 'velocity.h5')
    keep_open = False
    dset_names = 'timeseries'
    num_threads = 3
    block_shape = (256, 256)

    # extract one product to serve as a reference file
    ds = gdal.Translate(dolphin_ref_tif, unw_files[0])
    ds = None

    # get list of times WRT to the reference time
    with h5py.File(og_ts_file, 'r') as f:
        ts_date_list = f['date'][:]

    x_arr = [
        datetime.datetime.strptime(
            date.decode('utf-8'), '%Y%m%d'
        )
        for date in ts_date_list
    ]
    x_arr = datetime_to_float(x_arr)

    # initiate dolphin file object
    writer = BackgroundRasterWriter(dolphin_vel_file,
        like_filename=dolphin_ref_tif)
    s = HDF5StackReader.from_file_list(
            [og_ts_file], dset_names=dset_names, keep_open=keep_open
        )

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

    readers = [s[0,:,:]]
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
    start_date = date12_list[0].split('_')[0]
    end_date = date12_list[-1].split('_')[-1]
    meta["DATA_TYPE"] = 'float32'
    meta["DATE12"] =  start_date + '_' + end_date
    meta["FILE_PATH"] = og_ts_file
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

    # generate mask file from unw phase field
    if inps.zero_mask is True:
        msk_file = os.path.join(os.path.dirname(og_ts_file),
            'combined_msk.h5')
        iargs = [stack_file, 'unwrapPhase', '-o', msk_file, '--nonzero']
        generate_mask.main(iargs)

        # mask TS file, since reference_point adds offset back in masked field
        iargs = [og_ts_file, '--mask', msk_file]
        mask.main(iargs)
        # capture masked TS file
        ts_file = os.path.join(inps.out_dir, "timeseries_msk.h5")

        # mask velocity file
        iargs = [vel_file, '--mask', msk_file]
        mask.main(iargs)

    print("Done.")
    return


if __name__ == "__main__":
    main(sys.argv[1:])
