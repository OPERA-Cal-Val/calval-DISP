#!/usr/bin/env python

import os
import glob
import rasterio
from rasterio.enums import Resampling
import numpy as np
import h5py
import datetime

from mintpy.utils import writefile

import argparse

def parse_arguments():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(description="Process and save raster differences into an HDF5 file.")
    parser.add_argument(
        "--dir1", required=True, help="Path to the first directory containing the extracted TS epochs."
    )
    parser.add_argument(
        "--dir2", required=True, help="Path to the second directory containing the extracted TS epochs."
    )
    parser.add_argument(
        "--output_h5_path", default='diff_overlap.h5', help="Path to the output HDF5 file."
    )
    return parser.parse_args()


def match_files(dir1, dir2, wildcard="*.tif"):
    """
    Matches files by name between two directories based on a wildcard pattern.
    Returns a list of tuples containing matched file paths.
    """
    files1 = sorted(
        glob.glob(f'{dir1}/{wildcard}'),
        key=lambda x: datetime.datetime.strptime(
            x.split("timeseries-")[-1][:8], '%Y%m%d'
        )
    )
    files1 = {os.path.basename(f): f for f in files1}
    files2 = sorted(
        glob.glob(f'{dir2}/{wildcard}'),
        key=lambda x: datetime.datetime.strptime(
            x.split("timeseries-")[-1][:8], '%Y%m%d'
        )
    )
    files2 = {os.path.basename(f): f for f in files2}
    return [
        (files1[name], files2[name]) for name in files1 if name in files2
    ]


def subtract_rasters(raster1_path, raster2_path):
    """
    Subtracts two rasters over their common extent, resampling to the finest resolution.
    Returns the result data, transform, nodata value, and CRS.
    """
    with rasterio.open(raster1_path) as src1, rasterio.open(raster2_path) as src2:
        # Calculate common bounds
        bounds1, bounds2 = src1.bounds, src2.bounds
        common_bounds = (
            max(bounds1.left, bounds2.left),
            max(bounds1.bottom, bounds2.bottom),
            min(bounds1.right, bounds2.right),
            min(bounds1.top, bounds2.top),
        )

        # Determine resolution and dimensions
        resolution = min(src1.res[0], src2.res[0])
        width = int((common_bounds[2] - common_bounds[0]) / resolution)
        height = int((common_bounds[3] - common_bounds[1]) / resolution)
        transform = rasterio.transform.from_bounds(
            *common_bounds, width=width, height=height
        )

        # Read and resample data
        data1 = src1.read(
            1,
            out_shape=(height, width),
            resampling=Resampling.nearest,
            window=rasterio.windows.from_bounds(*common_bounds, src1.transform),
        )
        data1 = np.where(data1 == 0, np.nan, data1)
        data2 = src2.read(
            1,
            out_shape=(height, width),
            resampling=Resampling.nearest,
            window=rasterio.windows.from_bounds(*common_bounds, src2.transform),
        )
        data2 = np.where(data2 == 0, np.nan, data2)
        data2 = np.where(data1 == np.nan, np.nan, data2)

        # Subtract rasters
        result_data = data1 - data2

        # Handle nodata values
        nodata_value = src1.nodata or src2.nodata
        if nodata_value is not None:
            result_data[
                np.isclose(data1, nodata_value) | np.isclose(data2, nodata_value)
            ] = nodata_value

        # set nan to 0
        result_data = np.nan_to_num(result_data, nan=0)

        return result_data, transform, nodata_value, src1.crs


def save_to_hdf5(output_rasters, output_h5_path):
    """
    Saves a list of raster outputs to an HDF5 file with metadata.
    """
    # Get dimensions from the first raster
    first_raster = output_rasters[0][0]
    height, width = first_raster.shape
    num_bands = len(output_rasters)
    geotransform = output_rasters[0][1]
    gdal_transform = (geotransform[2], geotransform[0], geotransform[1],
                  geotransform[5], geotransform[3], geotransform[4])
    date_list = [i[4] for i in output_rasters]
    pbase = np.zeros(num_bands, dtype=np.float32)
    print('geotransform', geotransform)
    print('gdal_transform', gdal_transform)
    print('date_list', date_list)
    dates = np.array(date_list, dtype=np.string_)

    # initiate file
    ds_name_dict = {
        "date": [dates.dtype, (num_bands,), dates],
        "bperp": [np.float32, (num_bands,), pbase],
        "timeseries": [np.float32, (num_bands, height, width), None],
    }

    meta = {}
    meta["FILE_TYPE"] = "timeseries"
    meta["UNIT"] = "m"
    meta["LENGTH"] = height
    meta["WIDTH"] = width
    meta["X_FIRST"] = geotransform[2]
    meta["Y_FIRST"] = geotransform[5]
    meta["X_STEP"] = geotransform[0]
    meta["Y_STEP"] = geotransform[4]
    meta["X_UNIT"] = meta["Y_UNIT"] = "meters"

    writefile.layout_hdf5(output_h5_path, ds_name_dict, metadata=meta)

    with h5py.File(output_h5_path, "a") as h5f:
        for i, (data, transform, nodata_value, crs, name) in \
            enumerate(output_rasters):
            h5f["timeseries"][i] = data

        print(f"Saved all rasters to {output_h5_path}")


def process_and_save(dir1, dir2, output_h5_path, wildcard="*.tif"):
    """
    Matches rasters in two directories, subtracts them,
    and saves the results in an HDF5 file.
    """
    raster_pairs = match_files(dir1, dir2, wildcard=wildcard)
    output_rasters = []

    for raster1, raster2 in raster_pairs:
        result_data, transform, nodata_value, crs = subtract_rasters(raster1, raster2)
        date = os.path.basename(raster1).split("timeseries-")[-1][:8]
        nodata_value = 0
        output_rasters.append((result_data, transform, nodata_value, crs, date))

    save_to_hdf5(output_rasters, output_h5_path)


# Example Usage
if __name__ == "__main__":
    #dir1 = "/u/duvel-d0/ssangha/DISP_work/v9/dolphin/F11116/extracted_DISP"
    #dir2 = "/u/duvel-d0/ssangha/DISP_work/v9/dolphin/F11115/extracted_DISP"
    #output_h5_path = "/u/duvel-d0/ssangha/DISP_work/v9/dolphin/F11116-F11115/DISP/diff_overlap.h5"
    args = parse_arguments()

    # Use the arguments
    dir1 = args.dir1
    dir2 = args.dir2
    output_h5_path = args.output_h5_path

    process_and_save(dir1, dir2, output_h5_path)

