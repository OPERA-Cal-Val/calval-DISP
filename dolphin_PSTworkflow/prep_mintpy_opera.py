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
from tqdm import tqdm

# Related third-party imports
import affine
import h5py
import numpy as np
import pyproj
import rasterio
from rasterio import CRS
from rasterio.warp import reproject, Resampling

# Add the src directory to sys.path
sys.path.append(str(Path(__file__).parent / 'src'))

# Local application/library specific imports
from mintpy.cli import temporal_average, mask, generate_mask, \
    timeseries2velocity
from mintpy.reference_point import reference_point_attribute
from mintpy.utils import arg_utils, ptime, readfile, writefile
from mintpy.utils import utils as ut
from mintpy.utils.utils0 import calc_azimuth_from_east_north_obs
from mintpy.utils.utils0 import azimuth2heading_angle
from pst_dolphin_utils import create_external_files, get_raster_bounds, \
    get_raster_crs, full_suffix, get_raster_xysize, get_raster_gt, \
    load_gdal, warp_to_match
from tile_mate import get_raster_from_tiles
from tile_mate.stitcher import DATASET_SHORTNAMES
from opera_utils import get_dates
from dem_stitcher.stitcher import stitch_dem

OPERA_DATASET_ROOT = './'


####################################################################################
EXAMPLE = """example:

  prep_mintpy.py
      -m pst_output/static_CSLCs/
      -c "pst_output/dolphin_output/stitched_interferograms/*.zeroed.cor.tif"
      -u "pst_output/dolphin_output/stitched_interferograms/*.unw.zeroed.tif"
      --geom-dir pst_output/dolphin_output/stitched_interferograms/geometry
      --ref-lalo '19.2485991551617 -155.32285148610057'
      --single-reference
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
        "-b",
        "--baseline-dir",
        dest="baselineDir",
        type=str,
        default=None,
        help="baseline directory (default: %(default)s).",
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
        "--single-reference",
        dest="single_reference",
        action="store_true",
        help=(
            "Indicate that all the unwrapped ifgs are single reference, which allows"
            " use to create the timeseries.h5 file directly without inversion."
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
    num_file,
    file_list,
    ref_y,
    ref_x,
    water_mask,
    date12_list,
    phase2range=1,
    unw_file=False
    ):
    """Prepare h5 file for input stack of layers"""
    # initiate file
    writefile.layout_hdf5(fname, ds_name_dict, metadata=meta)

    # writing data to HDF5 file
    print("writing data to HDF5 file {} with a mode ...".format( \
          fname))
    with h5py.File(fname, "a") as f:
        prog_bar = ptime.progressBar(maxValue=num_file)
        for i,file_inc in enumerate(file_list):
            # read data using gdal
            # if unw file, cannot read as masked array
            # necessary to avoid errors with MintPy
            if unw_file is True:
                data = load_gdal(file_inc)
            else:
                data = load_gdal(file_inc, masked=True)
            # apply reference point
            data -= np.nan_to_num(data[ref_y, ref_x])
            # if unw file, convert nans to 0
            # necessary to avoid errors with MintPy velocity fitting
            if unw_file is True:
                data = np.nan_to_num(data)
            # also apply water mask
            data = data * water_mask
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
    baseline_dir=None,
):
    """Prepare the timeseries file."""
    print("-" * 50)
    print("preparing timeseries file: {}".format(outfile))

    # copy metadata to meta
    meta = {key: value for key, value in metadata.items()}

    # pass lyr to disp conversion factor, which depends on the product version
    if track_version == 0.3:
        disp_lyr_name = 'unwrapped_phase'
        phase2range = -1 * float(meta["WAVELENGTH"]) / (4.0 * np.pi)
    if track_version >= 0.4:
        disp_lyr_name = 'displacement'
        phase2range = 1

    # grab date list from the filename
    date12_list = _get_date_pairs(unw_files)
    num_file = len(unw_files)
    print("number of unwrapped interferograms: {}".format(num_file))

    date_pairs = [dl.split("_") for dl in date12_list]
    date_list = sorted(set(itertools.chain.from_iterable(date_pairs)))
    # ref_date = date12_list[0].split("_")[0]
    # date_list = [ref_date] + [date12.split("_")[1] for date12 in date12_list]
    num_date = len(date_list)
    print("number of acquisitions: {}\n{}".format(num_date, date_list))

    # baseline info
    pbase = np.zeros(num_date, dtype=np.float32)
    if baseline_dir is not None:
        raise NotImplementedError

    # size info
    cols, rows = get_raster_xysize(unw_files[0])

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

    # initiate file attributes
    meta["FILE_TYPE"] = "timeseries"
    meta["UNIT"] = "m"
    # meta["REF_DATE"] = ref_date # might not be the first date!
    ##################################################
    # pass reference point
    # if not specified, chose reference point from max value in coherence file
    coord = ut.coordinate(meta)
    if ref_lalo is not None:
        ref_lat, ref_lon = ref_lalo.split()
        # get value at coordinate
        (ref_y,
         ref_x) = coord.geo2radar(np.array(float(ref_lat)),
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

    # update metadata to reflect reference point
    ref_meta = reference_point_attribute(meta, y=ref_y, x=ref_x)
    meta['REF_LAT'] = ref_meta['REF_LAT']
    meta['REF_LON'] = ref_meta['REF_LON']
    meta['REF_Y'] = ref_meta['REF_Y']
    meta['REF_X'] = ref_meta['REF_X']

    # write TS containing displacement layers to file
    save_stack(outfile, ds_name_dict, meta, num_file, unw_files, ref_y, ref_x,
        water_mask, date12_list, phase2range, unw_file=True)

    # record output
    all_outputs = [outfile]

    # loop through and create files for each correction layer
    correction_layers = ['tropospheric_delay', 'ionospheric_delay', \
                         'solid_earth_tide', 'plate_motion']
    for lyr in correction_layers:
        # get correction layers
        corr_fname = os.path.join(os.path.dirname(outfile), f'{lyr}.h5')
        corr_files = \
            [i.replace(f'{disp_lyr_name}', f'/corrections/{lyr}') \
             for i in unw_files]

        # write layers to file
        save_stack(corr_fname, ds_name_dict, meta, num_file, corr_files,
            ref_y, ref_x, water_mask, date12_list, phase2range)

        # record outputs
        all_outputs.append(corr_fname)

    # save short wavelength layers to file
    if track_version >= 0.4:
        shortwvl_lyr = 'short_wavelength_displacement'
        # get correction layers
        shortwvl_fname = os.path.join(os.path.dirname(outfile),
            f'{shortwvl_lyr}.h5')
        shortwvl_files = \
            [i.replace(f'{disp_lyr_name}', f'{shortwvl_lyr}') \
             for i in unw_files]

        # write layers to file
        save_stack(shortwvl_fname, ds_name_dict, meta, num_file,
            shortwvl_files, ref_y, ref_x, water_mask, date12_list,
            phase2range)

        # record output
        all_outputs.append(shortwvl_fname)

    return all_outputs


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


def prepare_stack(
    outfile,
    unw_files,
    cor_files,
    track_version,
    metadata,
    water_mask_file=None,
    ref_lalo=None,
    # baseline_dir=None,
):
    """Prepare the input unw stack."""
    print("-" * 50)
    print("preparing ifgramStack file: {}".format(outfile))
    # copy metadata to meta
    meta = {key: value for key, value in metadata.items()}

    # get list of *.unw file
    num_pair = len(unw_files)
    unw_ext = full_suffix(unw_files[0])

    conv_factor = 1
    # >=v0.4 increment in units of m, must be converted to phs
    if track_version >= 0.4:
        conv_factor = -1 * (4.0 * np.pi) / float(metadata["WAVELENGTH"])

    print(unw_files)
    print(f"number of unwrapped interferograms: {num_pair}")
    print(f"number of correlation files: {len(cor_files)}")
    print(cor_files)

    # get list of *.unw.conncomp file
    cc_files = \
        [i.replace('unwrapped_phase', 'connected_component_labels') \
         for i in unw_files]

    print(f"number of connected components files: {len(cc_files)}")

    if len(cc_files) != len(unw_files) or len(cor_files) != len(unw_files):
        print(
            "the number of *.unw and *.unw.conncomp or *.cor files are NOT consistent"
        )
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
        "connectComponent": [
            np.float32,
            (num_pair, rows, cols),
            None,
        ],
    }

    # read water mask
    if water_mask_file is not None:
         water_mask = readfile.read(water_mask_file,
             datasetName='waterMask')[0]
    else:
        water_mask = np.ones((rows, cols), dtype=np.float32)

    # initiate HDF5 file
    meta["FILE_TYPE"] = "ifgramStack"
    writefile.layout_hdf5(outfile, ds_name_dict, metadata=meta)

    # writing data to HDF5 file
    print("writing data to HDF5 file {} with a mode ...".format(outfile))
    with h5py.File(outfile, "a") as f:
        prog_bar = ptime.progressBar(maxValue=num_pair)
        for i, (unw_file, cor_file, cc_file) in enumerate(
            zip(unw_files, cor_files, cc_files)
        ):
            # read/write *.unw file
            f["unwrapPhase"][i] = load_gdal(
                unw_file, masked=True) * water_mask * conv_factor

            # read/write *.cor file
            f["coherence"][i] = load_gdal(
                cor_file, masked=True) * water_mask

            # read/write *.unw.conncomp file
            f["connectComponent"][i] = load_gdal(
                cc_file, masked=True) * water_mask

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
                                f'{msk_file} to inform selection of new point.')

    # extract correction layers
    # Remove dict entries not relevant to correction layers
    for key in ['coherence', 'connectComponent']:
        ds_name_dict.pop(key, None)

    # loop through and create files for persistent scatterer
    ps_files = \
            [i.replace('unwrapped_phase', 'persistent_scatterer_mask') \
             for i in unw_files]
    ps_fname = os.path.join(os.path.dirname(outfile), 
        'persistent_scatterer_mask.h5')
    prepare_average_stack(ps_fname, ps_files, water_mask, 'mask', meta)

    # loop through and create files for temporal coherence
    tempcoh_files = \
            [i.replace('unwrapped_phase', 'temporal_coherence') \
             for i in unw_files]
    tempcoh_fname = os.path.join(os.path.dirname(outfile),
        'temporalCoherence.h5')
    prepare_average_stack(tempcoh_fname, tempcoh_files, water_mask,
                          'temporalCoherence', meta)

    return


def main(iargs=None):
    """Run the preparation functions."""
    inps = cmd_line_parse(iargs)

    unw_files = sorted(glob.glob(inps.unw_file_glob))

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
    cor_files = \
        [f'NETCDF:"{i}":interferometric_correlation' for i in unw_files]
    unw_files = \
        [f'NETCDF:"{i}":{disp_lyr_name}' for i in unw_files]
    print(f"Found {len(cor_files)} correlation files")

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
        cor_files=cor_files,
        track_version=track_version,
        metadata=meta,
        water_mask_file=inps.water_mask_file,
        ref_lalo=inps.ref_lalo,
        # baseline_dir=inps.baseline_dir,
    )

    # prepare TS file
    og_ts_file = os.path.join(inps.out_dir, "timeseries.h5")
    geom_file = os.path.join(inps.out_dir, "geometryGeo.h5")
    all_outputs = [og_ts_file]
    if inps.single_reference:
        # time-series (if inputs are all single-reference)
        all_outputs = prepare_timeseries(
            outfile=og_ts_file,
            unw_files=unw_files,
            track_version=track_version,
            metadata=meta,
            water_mask_file=inps.water_mask_file,
            ref_lalo=inps.ref_lalo,
            # baseline_dir=inps.baseline_dir,
        )

    mintpy_prepare_geometry(geom_file, geom_dir=inps.geom_dir, metadata=meta,
                            water_mask_file=inps.water_mask_file)

    # generate velocity fit
    vel_file = os.path.join(inps.out_dir, "velocity.h5")
    iargs = [og_ts_file, '-o', vel_file]
    timeseries2velocity.main(iargs)

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
