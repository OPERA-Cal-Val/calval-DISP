#!/usr/bin/env python3
"""PST miscellaneous utils for time-series management and tropo routines"""

import os
from datetime import datetime as dt, timedelta

import asf_search as asf
import h5py
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import rasterio
from affine import Affine
from eof.download import download_eofs
from osgeo import osr
from pst_dolphin_utils import load_gdal
from RAiDER.delay import tropo_delay
from RAiDER.llreader import BoundingBox
from RAiDER.losreader import Raytracing
from RAiDER.models.hrrr import HRRR
from RAiDER.processWM import prepareWeatherModel
from rasterio.warp import reproject, Resampling
from shapely import wkt

def calculate_cumulative_displacement(
    date, date_list, water_mask, mask_dict, lyr_name,
    rows, cols, ref_y, ref_x, G, phase2range,
    apply_tropo_correction, work_dir, median_height
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
        reflyr_name = file.split(':')[-1]
        file = file.replace(reflyr_name, lyr_name)
        print(f'{i} {ref_date} to {sec_date}: {file}')	

        # Read displacement data
        data = load_gdal(file, masked=True)
        data *= phase2range
        data = data.astype(np.float32)  # Convert to numeric type

        # Apply tropospheric correction if requested
        if apply_tropo_correction and work_dir:
            print(f"\nApplying tropospheric correction to pair {ref_date}-{sec_date}")
            
            # Extract original NetCDF filename from the GDAL path string
            # Format is 'NETCDF:"path/to/file.nc":displacement'
            nc_file = file.split('"')[1]
            print('nc file: ', nc_file)
            
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
                    if 'unwrapper_mask' in nc.keys():
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
                if 'unwrapper_mask' in locals():
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
                if 'unwrapper_mask' in locals():
                    del unwrapper_mask

                data -= calculated_tropo_delay      # unit: meter 

            except Exception as e:
                print(f"Warning: Tropospheric correction failed for {ref_date}_{sec_date}: {str(e)}")
                print("Continuing with uncorrected data...")

        # Apply reference point correction
        if ref_y is not None and ref_x is not None:
            data -= np.nan_to_num(data[ref_y, ref_x])
        
        # mask by specified dict of thresholds 
        for dict_key in mask_dict.keys():
            mask_lyr = file.replace(lyr_name, dict_key)
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
    # If RAiDER output is two-way delay
    diff_delays_ds = (delays[0] - delays[1])
    # If RAiDER output is one-way delay
    # diff_delays_ds = (delays[0] - delays[1]) * 2
    
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
    Search for Sentinel-1 scenes and return the sensor (S1A or S1B)
    of the first found scene.
    
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
