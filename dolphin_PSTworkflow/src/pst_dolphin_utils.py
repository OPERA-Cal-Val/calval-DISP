#!/usr/bin/env python3
"""Utils for accessing PST CSLC bursts and executing dolphin workflow"""

import logging
import os
from os import fspath
from pathlib import Path
from typing import (
    NamedTuple,
    Union
)

import affine
import numpy as np
import rasterio
from osgeo import gdal, gdal_array
from pyproj import CRS, Proj
from rasterio.warp import Resampling, reproject

from dem_stitcher.stitcher import stitch_dem
from tile_mate import get_raster_from_tiles
from tile_mate.stitcher import DATASET_SHORTNAMES

logger = logging.getLogger(__name__)


class Bbox(NamedTuple):
    """Bounding box named tuple, defining extent in cartesian coordinates.

    Usage:

        Bbox(left, bottom, right, top)

    Attributes
    ----------
    left : float
        Left coordinate (xmin)
    bottom : float
        Bottom coordinate (ymin)
    right : float
        Right coordinate (xmax)
    top : float
        Top coordinate (ymax)

    """

    left: float
    bottom: float
    right: float
    top: float


def create_external_files(lyr_name,
                          ref_file,
                          out_bounds,
                          crs,
                          out_dir='./',
                          maskfile=False,
                          demfile=False):
    """Generate mask and/or DEM file"""
    # get lat/lon bounds to generate mask and dem
    utm_proj = Proj(crs.to_wkt())
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
    resize_col, resize_row = get_raster_xysize(ref_file)
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


def gdal_to_numpy_type(gdal_type: Union[str, int]) -> np.dtype:
    """Convert gdal type to numpy type."""
    if isinstance(gdal_type, str):
        gdal_type = gdal.GetDataTypeByName(gdal_type)
    return np.dtype(gdal_array.GDALTypeCodeToNumericTypeCode(gdal_type))


def get_raster_crs(filename) -> CRS:
    """Get the CRS from a file.

    Parameters
    ----------
    filename
        Path to the file to load.

    Returns
    -------
    CRS
        CRS.

    """
    ds = gdal.Open(fspath(filename))

    return CRS.from_wkt(ds.GetProjection())


def get_raster_bounds(
    filename = None, ds = None
) -> Bbox:
    """Get the (left, bottom, right, top) bounds of the image."""
    if ds is None:
        if filename is None:
            msg = "Must provide either `filename` or `ds`"
            raise ValueError(msg)
        ds = gdal.Open(fspath(filename))

    gt = ds.GetGeoTransform()
    xsize, ysize = ds.RasterXSize, ds.RasterYSize

    left, top = _apply_gt(gt=gt, x=0, y=0)
    right, bottom = _apply_gt(gt=gt, x=xsize, y=ysize)

    return Bbox(left, bottom, right, top)


def _apply_gt(
    ds=None, filename=None, x=None, y=None, inverse=False, gt=None
) -> tuple[float, float]:
    """Read the (possibly inverse) geotransform, apply to the x/y coordinates."""
    if gt is None:
        if ds is None:
            ds = gdal.Open(fspath(filename))
            gt = ds.GetGeoTransform()
            ds = None
        else:
            gt = ds.GetGeoTransform()

    if inverse:
        gt = gdal.InvGeoTransform(gt)
    # Reference: https://gdal.org/tutorials/geotransforms_tut.html
    x = gt[0] + x * gt[1] + y * gt[2]
    y = gt[3] + x * gt[4] + y * gt[5]
    return x, y


def full_suffix(filename):
    """Get the full suffix of a filename, including multiple dots.

    Parameters
    ----------
    filename : str or Path
        path to file

    Returns
    -------
    str
        The full suffix, including multiple dots.

    Examples
    --------
        >>> full_suffix('test.tif')
        '.tif'
        >>> full_suffix('test.tar.gz')
        '.tar.gz'

    """
    fpath = Path(filename)
    return "".join(fpath.suffixes)


def get_raster_xysize(filename) -> tuple[int, int]:
    """Get the xsize/ysize of a GDAL-readable raster."""
    ds = gdal.Open(fspath(filename))
    xsize, ysize = ds.RasterXSize, ds.RasterYSize
    ds = None
    return xsize, ysize


def get_raster_gt(filename) -> list[float]:
    """Get the geotransform from a file.

    Parameters
    ----------
    filename
        Path to the file to load.

    Returns
    -------
    List[float]
        6 floats representing a GDAL Geotransform.

    """
    ds = gdal.Open(fspath(filename))
    return ds.GetGeoTransform()


def get_raster_driver(filename) -> str:
    """Get the GDAL driver `ShortName` from a file.

    Parameters
    ----------
    filename
        Path to the file to load.

    Returns
    -------
    str
        Driver name.

    """
    ds = gdal.Open(fspath(filename))
    return ds.GetDriver().ShortName


def get_raster_nodata(filename, band: int = 1):
    """Get the nodata value from a file.

    Parameters
    ----------
    filename
        Path to the file to load.
    band : int, optional
        Band to get nodata value for, by default 1.

    Returns
    -------
    Optional[float]
        Nodata value, or None if not found.

    """
    ds = gdal.Open(fspath(filename))
    return ds.GetRasterBand(band).GetNoDataValue()


def load_gdal(
    filename,
    *,
    band = None,
    subsample_factor = 1,
    overview = None,
    rows = None,
    cols = None,
    masked: bool = False,
) -> np.ndarray | np.ma.MaskedArray:
    """Load a gdal file into a numpy array.

    Parameters
    ----------
    filename : str or Path
        Path to the file to load.
    band : int, optional
        Band to load. If None, load all bands as 3D array.
    subsample_factor : int or tuple[int, int], optional
        Subsample the data by this factor. Default is 1 (no subsampling).
        Uses nearest neighbor resampling.
    overview: int, optional
        If passed, will load an overview of the file.
        Raster must have existing overviews, or ValueError is raised.
    rows : slice, optional
        Rows to load. Default is None (load all rows).
    cols : slice, optional
        Columns to load. Default is None (load all columns).
    masked : bool, optional
        If True, return a masked array using the raster's `nodata` value.
        Default is False.

    Returns
    -------
    arr : np.ndarray or np.ma.MaskedArray
        Array of shape (bands, y, x) or (y, x) if `band` is specified,
        where y = height // subsample_factor and x = width // subsample_factor.

    """
    ds = gdal.Open(fspath(filename))
    nrows, ncols = ds.RasterYSize, ds.RasterXSize

    if overview is not None:
        # We can handle the overviews most easily
        bnd = ds.GetRasterBand(band or 1)
        ovr_count = bnd.GetOverviewCount()
        if ovr_count > 0:
            idx = ovr_count + overview if overview < 0 else overview
            out = bnd.GetOverview(idx).ReadAsArray()
            bnd = ds = None
            return out
        logger.warning(f"Requested {overview = }, but none found for {filename}")

    # if rows or cols are not specified, load all rows/cols
    rows = slice(0, nrows) if rows in (None, slice(None)) else rows
    cols = slice(0, ncols) if cols in (None, slice(None)) else cols
    # Help out mypy:
    assert rows is not None
    assert cols is not None

    dt = gdal_to_numpy_type(ds.GetRasterBand(1).DataType)

    if isinstance(subsample_factor, int):
        subsample_factor = (subsample_factor, subsample_factor)

    xoff, yoff = int(cols.start), int(rows.start)
    row_stop = min(rows.stop, nrows)
    col_stop = min(cols.stop, ncols)
    xsize, ysize = int(col_stop - cols.start), int(row_stop - rows.start)
    if xsize <= 0 or ysize <= 0:
        msg = (
            f"Invalid row/col slices: {rows}, {cols} for file {filename} of size"
            f" {nrows}x{ncols}"
        )
        raise IndexError(msg)
    nrows_out, ncols_out = (
        ysize // subsample_factor[0],
        xsize // subsample_factor[1],
    )

    # Read the data, and decimate if specified
    resamp = gdal.GRA_NearestNeighbour
    if band is None:
        count = ds.RasterCount
        out = np.empty((count, nrows_out, ncols_out), dtype=dt)
        ds.ReadAsArray(xoff, yoff, xsize, ysize, buf_obj=out, resample_alg=resamp)
        if count == 1:
            out = out[0]
    else:
        out = np.empty((nrows_out, ncols_out), dtype=dt)
        bnd = ds.GetRasterBand(band)
        bnd.ReadAsArray(xoff, yoff, xsize, ysize, buf_obj=out, resample_alg=resamp)

    if not masked:
        return out
    # Get the nodata value
    nd = get_raster_nodata(filename)
    if nd is not None and np.isnan(nd):
        return np.ma.masked_invalid(out)
    else:
        return np.ma.masked_equal(out, nd)


def warp_to_match(
    input_file,
    match_file,
    output_file = None,
    resample_alg: str = "near",
    output_format = None,
) -> Path:
    """Reproject `input_file` to align with the `match_file`.

    Uses the bounds, resolution, and CRS of `match_file`.

    Parameters
    ----------
    input_file
        Path to the image to be reprojected.
    match_file
        Path to the input image to serve as a reference for the reprojected image.
        Uses the bounds, resolution, and CRS of this image.
    output_file
        Path to the output, reprojected image.
        If None, creates an in-memory warped VRT using the `/vsimem/` protocol.
    resample_alg: str, optional, default = "near"
        Resampling algorithm to be used during reprojection.
        See https://gdal.org/programs/gdalwarp.html#cmdoption-gdalwarp-r for choices.
    output_format: str, optional, default = None
        Output format to be used for the output image.
        If None, gdal will try to infer the format from the output file extension, or
        (if the extension of `output_file` matches `input_file`) use the input driver.

    Returns
    -------
    Path
        Path to the output image.
        Same as `output_file` if provided, otherwise a path to the in-memory VRT.

    """
    bounds = get_raster_bounds(match_file)
    crs_wkt = get_raster_crs(match_file).to_wkt()
    gt = get_raster_gt(match_file)
    resolution = (gt[1], gt[5])

    if output_file is None:
        output_file = f"/vsimem/warped_{Path(input_file).stem}.vrt"
        logger.debug(f"Creating in-memory warped VRT: {output_file}")

    if output_format is None and Path(input_file).suffix == Path(output_file).suffix:
        output_format = get_raster_driver(input_file)

    options = gdal.WarpOptions(
        dstSRS=crs_wkt,
        format=output_format,
        xRes=resolution[0],
        yRes=resolution[1],
        outputBounds=bounds,
        outputBoundsSRS=crs_wkt,
        resampleAlg=resample_alg,
    )
    gdal.Warp(
        fspath(output_file),
        fspath(input_file),
        options=options,
    )

    return Path(output_file)
