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
from osgeo import gdal
import rasterio
from tqdm import tqdm
import networkx as nx

# Add the src directory to sys.path
sys.path.append(str(Path(__file__).resolve().parents[1] / 'src'))

# Local application/library-specific imports
from mintpy.cli import (
    generate_mask,
    mask,
    temporal_average
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
from pst_ts_utils import calculate_cumulative_displacement
from tile_mate.stitcher import DATASET_SHORTNAMES

OPERA_DATASET_ROOT = './'


####################################################################################
EXAMPLE = """example:

  prep_mintpy.py
      -m pst_output/static_CSLCs/
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
        "-g",
        "--geom-dir",
        default="./geometry",
        help="Geometry directory (default: %(default)s).\n",
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
        "--tropo-correction",
        dest="tropo_correction",
        action="store_true",
        help="Apply tropospheric correction using HRRR weather model"
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
    reflyr_name = full_suffix(file_list[0])
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
    metadata,
    water_mask_file=None,
    ref_lalo=None,
    apply_tropo_correction=False,
    median_height=50,
    work_dir=None,
):
    """Prepare the timeseries file."""
    print("-" * 50)
    print("preparing timeseries file: {}".format(outfile))

    # get file ext
    unw_ext = full_suffix(unw_files[0])

    # copy metadata to meta
    meta = {key: value for key, value in metadata.items()}
    phase2range = 1

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
    if water_mask_file:
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
    meta['REF_LAT'] = ref_meta['REF_LAT']
    meta['REF_LON'] = ref_meta['REF_LON']
    meta['REF_Y'] = ref_meta['REF_Y']
    meta['REF_X'] = ref_meta['REF_X']

    # set dictionary that will be used to mask TS by specified thresholds
    mask_dict = {}

    print(f"Writing data to HDF5 file {outfile}")
    writefile.layout_hdf5(outfile, ds_name_dict, metadata=meta)

    reflyr_name = '.unw.tif'
    with h5py.File(outfile, "a") as f:
        prog_bar = ptime.progressBar(maxValue=num_date)
        # First date is always zero
        f["timeseries"][0] = np.zeros((rows, cols), dtype=np.float32)
        # Calculate cumulative displacement for each date
        for i, date in enumerate(date_list[1:], start=1):
            displacement = calculate_cumulative_displacement(
                date, date_list, water_mask, mask_dict, reflyr_name,
                rows, cols, ref_y, ref_x, G, phase2range,
                apply_tropo_correction, work_dir, median_height)
            if displacement is not None:
                # writing timeseries to the date
                f["timeseries"][i] = displacement
            prog_bar.update(i, suffix=date)
        
        prog_bar.close()

    all_outputs = [outfile]

    return all_outputs, ref_meta


def mintpy_prepare_geometry(outfile, geom_dir, metadata,
                            water_mask_file=None):
    """Prepare the geometry file."""
    print("-" * 50)
    print(f"preparing geometry file: {outfile}")

    geom_path = Path(geom_dir)
    # copy metadata to meta
    meta = {key: value for key, value in metadata.items()}
    meta["FILE_TYPE"] = "geometry"

    file_to_path = {
        "los_east": geom_path / "los_east.tif",
        "los_north": geom_path / "los_north.tif",
        "height": geom_path / "height.tif",
        "shadowMask": geom_path / "layover_shadow_mask.tif",
    }

    if water_mask_file:
        file_to_path["waterMask"] = water_mask_file

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
            print(f"Skipping {fname}: {e}")

    # Compute the azimuth and incidence angles from east/north coefficients
    azimuth_angle, east, north = get_azimuth_ang(dsDict)
    dsDict["azimuthAngle"] = azimuth_angle

    up = np.sqrt(1 - east**2 - north**2)
    incidence_angle = np.rad2deg(np.arccos(up))
    dsDict["incidenceAngle"] = incidence_angle

    writefile.write(dsDict, outfile, metadata=meta)
    return outfile


def prepare_average_stack(outfile, infiles, water_mask, file_type, metadata):
    """Average and export specified layers"""
    print("-" * 50)
    print("preparing average of stack: {}".format(outfile))

    # copy metadata to meta
    meta = {key: value for key, value in metadata.items()}
    meta["FILE_TYPE"] = file_type # "mask" "temporalCoherence"
    meta["UNIT"] = "1"

    # read data using gdal
    data = []
    for infile in infiles:
        data.append(load_gdal(infile, masked=True) * water_mask)
    data = np.nanmean(data, axis=0)

    # write to HDF5 file
    writefile.write(data, outfile, metadata=meta)
    return outfile


def main(iargs=None):
    """Run the preparation functions."""
    inps = cmd_line_parse(iargs)

    unw_files = sorted(
        glob.glob(inps.unw_file_glob),
        key=lambda x: datetime.datetime.strptime(
            x.split('_')[-1][:8], '%Y%m%d'
        )
    )
    print(f"Found {len(unw_files)} input tif files")
    print('unw_files', unw_files)

    # filter input by specified dates
    if inps.startDate is not None:
        startDate = int(inps.startDate)
        filtered_unw_files = []
        for i in enumerate(unw_files):
            sec_date = int(os.path.basename(i[1]).split('_')[1][:8])
            if sec_date >= startDate:
                filtered_unw_files.append(i[1])
        unw_files = filtered_unw_files

    if inps.endDate is not None:
        endDate = int(inps.endDate)
        filtered_unw_files = []
        for i in unw_files:
            sec_date = int(os.path.basename(i).split('_')[1][:8])
            if sec_date <= endDate:
                filtered_unw_files.append(i)
        unw_files = filtered_unw_files
    
    date12_list = _get_date_pairs(unw_files)
    print(f"Found {len(unw_files)} unwrapped files")

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
    # check if mask file already generated through dolphin
    else:
        msk_file = geometry_dir.glob('*_mask.tif')
        msk_file = [i for i in msk_file]
        if len(msk_file) > 0:
            inps.water_mask_file = msk_file[0]
            print(f"Found water mask file {inps.water_mask_file}")

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

    # get median height for tropo estimate
    with rasterio.open(height_file) as src:
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

    # prepare TS file
    og_ts_file = os.path.join(inps.out_dir, "timeseries.h5")
    geom_file = os.path.join(inps.out_dir, "geometryGeo.h5")
    all_outputs = [og_ts_file]
    ref_meta = None
    all_outputs, ref_meta = prepare_timeseries(
        outfile=og_ts_file,
        unw_files=unw_files,
        metadata=meta,
        water_mask_file=inps.water_mask_file,
        ref_lalo=inps.ref_lalo,
        apply_tropo_correction=inps.tropo_correction,
        median_height=median_height,
        work_dir=inps.work_dir if inps.tropo_correction else None,
    )

    mintpy_prepare_geometry(geom_file, geom_dir=inps.geom_dir, metadata=meta,
                            water_mask_file=inps.water_mask_file)

    # write dolphin velocity fit to h5
    dolphin_vel_file = os.path.join(os.path.dirname(unw_files[0]),
        'velocity.tif')
    vel_file = os.path.join(inps.out_dir, 'velocity.h5')

    # update metadata field
    meta["DATA_TYPE"] = 'float32'
    meta["DATE12"] = date12_list[-1]
    meta["FILE_PATH"] = og_ts_file
    meta["FILE_TYPE"] = 'velocity'
    meta["NO_DATA_VALUE"] = 'none'
    meta["PROCESSOR"] = 'dolphin'
    meta["REF_DATE"] = date12_list[-1].split('_')[0]
    meta["START_DATE"] = date12_list[-1].split('_')[0]
    meta["END_DATE"] = date12_list[-1].split('_')[-1]
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

    print("Done.")
    return


if __name__ == "__main__":
    main(sys.argv[1:])
