#!/usr/bin/env python3
"""Utils for accessing PST CSLC bursts and executing dolphin workflow"""

# Standard library imports
import os

# Related third-party imports
import affine
import pyproj
import rasterio
from rasterio.warp import reproject, Resampling

# Local application/library specific imports
from dolphin import io
from tile_mate import get_raster_from_tiles
from tile_mate.stitcher import DATASET_SHORTNAMES
from dem_stitcher.stitcher import stitch_dem


def create_external_files(lyr_name,
                          ref_file,
                          out_bounds,
                          crs,
                          out_dir='./',
                          maskfile=False,
                          demfile=False):
    """Generate mask and/or DEM file"""
    # get lat/lon bounds to generate mask and dem
    utm_proj = pyproj.Proj(crs.to_wkt())
    e_lon, n_lat = utm_proj(out_bounds.left,
                            out_bounds.top,
                            inverse=True)
    w_lon, s_lat = utm_proj(out_bounds.right,
                            out_bounds.bottom,
                            inverse=True)
    geo_bounds = [e_lon, s_lat, w_lon, n_lat]
    # Get the affine transformation matrix
    with rasterio.open(ref_file) as src:
        reference_gt = src.transform
    resize_col, resize_row = io.get_raster_xysize(ref_file)
    #
    # download mask
    if maskfile:
        dat_arr, dat_prof = get_raster_from_tiles(geo_bounds,
                            tile_shortname=lyr_name)
        # fill permanent water body
        if lyr_name == 'esa_world_cover_2021':
            dat_arr[dat_arr == 80] = 0
            dat_arr[dat_arr != 0] = 1
        dat_arr = dat_arr.astype('byte')
        output_name = os.path.join(out_dir,
                      f'{lyr_name}_mask.tif')
        f_dtype = 'uint8'
        resampling_mode = Resampling.nearest
        print(f'mask file from source {lyr_name}\n')
    # download DEM
    if demfile:
        dst_area_or_point = 'Point'
        dst_ellipsoidal_height = True
        dat_arr, dat_prof = stitch_dem(geo_bounds,
                  dem_name=lyr_name,
                  dst_ellipsoidal_height=dst_ellipsoidal_height,
                  dst_area_or_point=dst_area_or_point)
        output_name = os.path.join(out_dir,
                      f'{lyr_name}_DEM.tif')
        f_dtype = 'float32'
        resampling_mode = Resampling.bilinear
        print(f'DEM file from source {lyr_name}\n')
    # resample
    with rasterio.open(output_name, 'w', 
                       height=resize_row, width=resize_col,
                       count=1,
                       dtype=f_dtype,
                       crs=crs,
                       transform=affine.Affine(*reference_gt)) as dst:
        reproject(source=dat_arr,
                  destination=rasterio.band(dst, 1),
                  src_transform=dat_prof['transform'],
                  src_crs=dat_prof['crs'],
                  dst_transform=reference_gt,
                  dst_crs=crs,
                  resampling=resampling_mode
                  )
    print(f'downloaded here: {output_name}\n')

    return output_name
