#!/usr/bin/env python3
"""Tools for accessing PST CSLC bursts and executing dolphin workflow"""

# Standard library imports
import argparse
import datetime
import subprocess
import sys
import time
from pathlib import Path

# Related third-party imports
import asf_search as asf
import boto3
import fsspec
import geopandas as gpd
import pandas as pd
from botocore import UNSIGNED
from botocore.config import Config
from shapely.geometry import box
from tqdm import tqdm

import requests
from io import BytesIO
import zipfile
import json

# Add the src directory to sys.path
sys.path.append(str(Path(__file__).parent / 'src'))

# Local application/library specific imports
from dolphin.io import get_raster_bounds, get_raster_crs
from dolphin.utils import prepare_geometry
from dolphin.workflows import _cli_config as dconfig
from dolphin.workflows import _cli_run as drun
from dolphin.workflows import CallFunc, config, stitching_bursts, unwrapping
from dolphin.workflows.displacement import timeseries
from dem_stitcher.stitcher import stitch_dem
from pst_dolphin_utils import create_external_files
from tile_mate import get_raster_from_tiles
from tile_mate.stitcher import DATASET_SHORTNAMES


def create_parser():
    """
        Run dolphin workflow with specified input parameters
    """
    parser = argparse.ArgumentParser(description='Access CSLCs '
                                     'from S3 bucket and run dolphin')
    parser.add_argument('-s', '--startdate', dest='start_date', type=str,
                        default=None, help='start date YYYYMMDD')
    parser.add_argument('-e', '--enddate', dest='end_date', type=str,
                        default=None, help='end date YYYYMMDD')
    parser.add_argument('-fr', '--fixedorrange', action='store_true',
                        dest='fixed_or_range',
                        help='If true, input dates treated as fixed dates '
                             'for pair-wise processing. '
                             'By default input dates are '
                             'treated as a range.')
    parser.add_argument('-fi', '--frameid', dest='frame_id', type=int,
                        default=None, help='Specify OPERA frame ID')
    parser.add_argument('-ti', '--trackid', dest='track_id', type=int,
                        default=None, help='Specify Sentinel track ID')
    parser.add_argument('-op', '--orbitpass', dest='orbit_pass', type=str,
                        default=None,
                        help='Specify orbit pass direction. Either: '
                             '1. Ascending, or '
                             '2. Descending')
    parser.add_argument('-ao', '--areaofinterest',
                        dest='area_of_interest', type=str,
                        default=None,
                        help='Specify area of interest as: '
                             '1. path to a valid shapefile, or '
                             '2. S N W E coordinates, ie ISCE convention')
    parser.add_argument('--strides', dest='strides',
                        type=str,
                        default='6 3',
                        help='Ifg option: Specify the (x y) '
                             'strides (decimation factor')
    parser.add_argument('--threadsperworker', dest='threads_per_worker',
                        type=int,
                        default=8,
                        help='Ifg option: specify num of threads per worker')
    parser.add_argument('--unwmethod', dest='unwrap_method',
                        type=str,
                        default='snaphu',
                        help='Unw option: specify an unw method: '
                             '1. snaphu '
                             '2. icu '
                             '3. phass ')
    parser.add_argument('--nparalleljobs', dest='n_parallel_jobs',
                        type=int,
                        default=2,
                        help='Unw option: '
                             'specify num of IFGs to unw in parallel')
    parser.add_argument('--nparalleltiles', dest='n_parallel_tiles',
                        type=int,
                        default=2,
                        help='Unw option: '
                             'specify num of tiles to unw in parallel')
    parser.add_argument('--ntiles', dest='ntiles', type=int,
                        default=4,
                        help='Unw option: specify num of tiles. '
                             'Suggested range is 4-16')
    parser.add_argument('-o', '--outdir', dest='out_dir', type=str,
                        default='./', help='Specify output directory')
    parser.add_argument('-verbose', '--verbose', action='store_true',
                        dest='verbose', help='Toggle verbose mode on')
    parser.add_argument('--water-mask-file', dest='water_mask_file', type=str,
                        default='esa_world_cover_2021',
                        help='Specify either path to valid water mask, or '
                             'download using one of the following '
                             f'data sources: {DATASET_SHORTNAMES}')
    parser.add_argument('--dem-file', dest='dem_file', type=str,
                        default='glo_30',
                        help='Specify either path to valid DEM, or '
                             'download using one of the following data '
                             'sources: srtm_v3, nasadem, glo_30, glo_90 '
                             'glo_90_missing')
    parser.add_argument('--correlation-threshold',
                        dest='correlation_threshold', type=float,
                        default=0.2,
                        help='Specify the correlation threshold '
                             'for TS mode')

    return parser


def cmd_line_parse(iargs=None):
    """
        Parser wrapper
    """
    parser = create_parser()
    user_inp = parser.parse_args(args=iargs)

    return user_inp


def filter_consistent_dates(df):
    """
        Filter out dates not common to each date
        Filter out bursts with less than 3 dates
    """
    # check if each burst has at least 3 dates
    list_of_bursts = list(set(df['operaBurstID'].to_list()))
    reject_bursts = []
    for i in list_of_bursts:
        all_dates = df['date'][df['operaBurstID'] == i].to_list()
        if len(all_dates) < 3:
            reject_bursts.append(i)
    df = df[~df['operaBurstID'].isin(reject_bursts)]

    # check if dates are common to all bursts
    list_of_dates = list(set(df['date'].to_list()))
    list_of_bursts = list(set(df['operaBurstID'].to_list()))
    reject_dates = []
    for i in list_of_dates:
        valid_burstlist = df['operaBurstID'][df['date'] == i].to_list()
        if len(valid_burstlist) != len(list_of_bursts):
            reject_dates.append(i)
    df = df[~df['date'].isin(reject_dates)]

    # one more time check if each burst has at least 2 dates
    list_of_bursts = list(set(df['operaBurstID'].to_list()))
    reject_bursts = []
    for i in list_of_bursts:
        all_dates = df['date'][df['operaBurstID'] == i].to_list()
        if len(all_dates) < 2:
            reject_bursts.append(i)
    df = df[~df['operaBurstID'].isin(reject_bursts)]

    return df


def create_inpt_txtfile(output_path, target_prefix, output_txtfile):
    """
        Record used PST CSLCs
    """
    # Get a list of files with the specified extension in the directory
    query_outputs = []
    for h5_file in output_path.iterdir():
        if h5_file.is_file() and target_prefix in str(h5_file):
            query_outputs.append(str(h5_file))

    # Write the list of files to the text file
    with output_txtfile.open(mode='w') as txtfile:
        txtfile.write('\n'.join(query_outputs))

    return


def map_bursts_to_frameids(repo_url, repo_zip, json_file):
    """
        Access json matching bursts to frame IDs
    """
    # URL of the ZIP file containing the JSON file
    repo_zip_url = repo_url + repo_zip

    # Download the ZIP file
    response = requests.get(repo_zip_url)
    zip_data = BytesIO(response.content)

    # Extract the JSON file from the ZIP archive
    with zipfile.ZipFile(zip_data, 'r') as zip_ref:
        # Assuming your JSON file is named 'data.json' within the ZIP
        json_data = zip_ref.read(json_file)

    # Load the JSON data
    burst_to_frame_json = json.loads(json_data.decode('utf-8'))['data']

    return burst_to_frame_json


def access_cslcs(inps=None):
    """
        Custom Dolphin workflow to access and leverage PST CSLCs
    """
    # Initialize variables and workspace
    csv_path = (
    's3://opera-provisional-products/DISP/DISP-S1/validation_data/'
    'getallbursts_table_validation_bursts_target_v1.0.csv'
    )
    start_date = inps.start_date
    end_date = inps.end_date
    fixed_or_range = inps.fixed_or_range
    frameid_number = inps.frame_id
    trackid_number = inps.track_id
    orbit_pass = inps.orbit_pass
    area_of_interest = inps.area_of_interest
    strides = inps.strides.split()
    threads_per_worker = inps.threads_per_worker
    n_parallel_jobs = inps.n_parallel_jobs
    n_parallel_tiles = inps.n_parallel_tiles
    ntiles = (inps.ntiles, inps.ntiles)
    unwrap_method = inps.unwrap_method
    water_mask_file = inps.water_mask_file
    dem_file = inps.dem_file
    correlation_threshold = inps.correlation_threshold
    if orbit_pass:
        if orbit_pass[:3].lower() == 'asc':
            orbit_pass = 'ASCENDING'
        if orbit_pass[:3].lower() == 'des':
            orbit_pass = 'DESCENDING'

    # initialize frameid and/or area of interest
    if not frameid_number and not area_of_interest:
        raise ValueError('Must specify valid input for at least either: '
                         '--fr (--frameid) or --ao (--areaofinterest)')
    if not frameid_number and not orbit_pass:
        raise ValueError('Must specify valid input for at least either: '
                         '--fr (--frameid) or --op (--orbitpass)')

    if area_of_interest:
        if Path(area_of_interest).is_file():
            area_of_interest = gpd.read_file(area_of_interest)
            area_of_interest = area_of_interest['geometry'][0]
        else:
            area_of_interest = area_of_interest.split()
            area_of_interest = [float(i) for i in area_of_interest]
            area_of_interest = box(area_of_interest[2],
                                   area_of_interest[0],
                                   area_of_interest[3],
                                   area_of_interest[1])  # (W,S,E,N)

    # set output dirs, and create them if necessary
    out_dir = Path(inps.out_dir).resolve()
    static_dir = out_dir.joinpath('static_CSLCs')
    cslc_dir = out_dir.joinpath('CSLCs')
    dolphin_dir = out_dir.joinpath('dolphin_output')
    stitched_ifg_path = dolphin_dir.joinpath('stitched_interferograms')
    for i in [out_dir, static_dir, cslc_dir, dolphin_dir, stitched_ifg_path]:
        if not i.exists():
            i.mkdir(parents=True, exist_ok=True)

    # Initiate dictionary to query for CSLCs
    asf_search_kwargs = {
        'dataset': 'OPERA-S1', 'processingLevel': 'CSLC',
        'intersectsWith': area_of_interest.wkt, 'flightDirection': orbit_pass}

    # convert the string to a datetime object
    date_format = "%Y%m%d"
    if start_date is not None:
        start_datetime = datetime.datetime.strptime(start_date, date_format)
        asf_search_kwargs['start'] = start_datetime
    if end_date is not None:
        end_datetime = datetime.datetime.strptime(end_date, date_format)
        end_datetime = end_datetime.replace(hour=23, minute=59, second=59)
        asf_search_kwargs['end'] = end_datetime

    # pass track id
    if trackid_number is not None:
        asf_search_kwargs['relativeOrbit'] = trackid_number

    # pass frame id
    if frameid_number is not None:
        asf_search_kwargs['asfFrame'] = frameid_number

    ## Return asf results to a dataframe
    results = asf.search(**asf_search_kwargs)
    print(f"Length of Results: {len(results)}")
    df = gpd.GeoDataFrame.from_features(results.geojson())

    # parse datetime
    df['date'] = pd.to_datetime(df['startTime'],
        format='%Y-%m-%dT%H:%M:%SZ').dt.strftime(date_format)

    # If specified, only capture input dates which match input
    # as opposed to treating them as a range
    if start_date is not None and end_date is not None:
        if fixed_or_range:
            # Filter dataframe by specified date range
            df = df[(df['date'] == start_date)
                    | (df['date'] == end_date)]

    # Get map of burst ID to frame ID
    repo_url = 'https://github.com/opera-adt/burst_db/releases/download/' + \
        'v0.3.1/'
    repo_zip = 'opera-s1-disp-burst-to-frame-0.3.1.json.zip'
    json_file = 'opera-s1-disp-burst-to-frame.json'
    burst_to_frame_json = map_bursts_to_frameids(repo_url, repo_zip,
        json_file)
    # get frame IDs
    df['frame_ids'] = df['operaBurstID'].str.lower().map(burst_to_frame_json)
    df['frame_ids'] = df['frame_ids'].apply(lambda x: x['frame_ids'])

    # determine and pass most common path if AOI specified, but no frame ID
    if area_of_interest and not frameid_number:
        all_frameids = df['frame_ids'].to_list()
        all_frameids = [item for sublist in all_frameids for item in sublist]
        # Use a dictionary to count occurrences
        count_dict = {}
        for item in all_frameids:
            count_dict[item] = count_dict.get(item, 0) + 1
        # Find the most common value
        frameid_number = max(count_dict, key=count_dict.get)
        print('Valid --fr (--frameid) not specified, but based on defined '
              f'--ao (--areaofinterest) of {area_of_interest}, the most '
              f'common frameid is {frameid_number}. The dataframe will be '
              'filtered by this value.')

    # filter by track
    if frameid_number is not None:
        df = df[df['frame_ids'].apply(lambda x: frameid_number in x)]

    # filter by only dates common to each burst
    # and remove bursts with less than 3 dates
    df = filter_consistent_dates(df)

    # Reset the index of the filtered DataFrame
    df = df.reset_index(drop=True)

    # download static layer just for the first dates
    list_of_bursts = list(set(df['operaBurstID'].to_list()))
    # search and download CLSC Static Layer files
    results_static = asf.search(
        operaBurstID=list_of_bursts,
        processingLevel='CSLC-STATIC')
    results_static.download(path=static_dir, processes=5)

    #
    # Produce burstwise IFGs over filtered subset
    # run dolphin by burst
    list_of_bursts = list(set(df['operaBurstID'].to_list()))
    for burst in tqdm(list_of_bursts, desc="Processing burst:"):
        print(f"{burst}")
        # set output path for burst
        dolphin_dir_burst = dolphin_dir.joinpath(f'{burst}')
        if not dolphin_dir_burst.exists():
            dolphin_dir_burst.mkdir(parents=True, exist_ok=True)
        # filter by burst
        filtered_df = df[df['operaBurstID'] == burst]
        filtered_df = filtered_df.reset_index(drop=True)
        cslc_txt = out_dir.joinpath(f'CSLC_list_{burst}.txt')
        static_txt = out_dir.joinpath(f'static_list_{burst}.txt')
        for index, row in tqdm(filtered_df.iterrows(),
                               desc="Downloading scene"):
            print(f"{row['date']}")
            # download CSLC
            url = row['url']
            local_filename = url.split("/")[-1]
            local_filename = cslc_dir.joinpath(local_filename)
            asf.download_url(url=url, path=cslc_dir)
            # record downloads in respective text-files
            target_prefix = 'OPERA_L2_CSLC-S1_T'
            create_inpt_txtfile(cslc_dir, target_prefix, cslc_txt)
            # capture static CSLCs
            if index == 0:
                target_prefix = 'OPERA_L2_CSLC-S1-STATIC_T'
                create_inpt_txtfile(static_dir, target_prefix, static_txt)

        # create config file and run dolphin
        yml_file = out_dir.joinpath(f'dolphin_config_{burst}.yaml')
        dconfig.create_config(outfile=str(yml_file),
                              slc_files=cslc_txt,
                              work_directory=str(dolphin_dir_burst),
                              subdataset='/data/VV',
                              strides=strides,
                              n_parallel_bursts=1,
                              threads_per_worker=threads_per_worker,
                              enable_gpu=True,
                              no_unwrap=True,
                              no_inversion=True,
                              zero_where_masked=True,
                              output_bounds=None)
        # run dolphin
        drun.run(str(yml_file))
        # delete CSLCs to save space
        process_cmd = f'rm -rf {str(cslc_dir)}/OPERA_L2_CSLC-S1_T*'
        subprocess.run(process_cmd, shell=True,
                       cwd=out_dir, capture_output=True,
                       text=True, check=True)
        # delete burst unw IFGs to save space
        process_cmd = f'rm -rf {str(dolphin_dir_burst)}/unwrapped'
        subprocess.run(process_cmd, shell=True,
                       cwd=out_dir, capture_output=True,
                       text=True, check=True)
        # delete burst unw IFGs to save space
        process_cmd = f'rm -rf {str(dolphin_dir_burst)}/timeseries'
        subprocess.run(process_cmd, shell=True,
                       cwd=out_dir, capture_output=True,
                       text=True, check=True)

    #
    # Stitch all ifgs by burst
    # set static paths for dolphin
    ifg_paths = list(dolphin_dir.glob('t*/interferograms/*.vrt'))
    coh_paths = list(dolphin_dir.glob('t*/linked_phase/' +
                     'temporal_coherence_average_*.tif'))
    shp_count_file_list = list(dolphin_dir.glob('t*/interferograms/' +
                               'shp_counts.tif'))
    # only access looked file if applicable
    if strides == ['1', '1']:
        ps_file_list = list(dolphin_dir.glob('t*/PS/ps_pixels.tif'))
        amp_dispersion_list = \
            list(dolphin_dir.glob('t*/PS/amp_dispersion.tif'))
    else:
        ps_file_list = list(dolphin_dir.glob('t*/PS/ps_pixels_looked.tif'))
        amp_dispersion_list = \
            list(dolphin_dir.glob('t*/PS/amp_dispersion_looked.tif'))
    cfg_obj = config.DisplacementWorkflow.from_yaml(yml_file)
    (
        stitched_ifg_paths,
        stitched_cor_paths,
        stitched_temp_coh_file,
        stitched_ps_file,
        stitched_amp_dispersion_file,
        stitched_shp_count_file,
    ) = stitching_bursts.run(
        ifg_file_list=ifg_paths,
        temp_coh_file_list=coh_paths,
        ps_file_list=ps_file_list,
        amp_dispersion_list=amp_dispersion_list,
        shp_count_file_list=shp_count_file_list,
        stitched_ifg_dir=stitched_ifg_path,
        output_options=cfg_obj.output_options,
        file_date_fmt='%Y%m%d',
        corr_window_size=(11, 11),
    )
    #
    # generate mask file
    geometry_dir = stitched_ifg_path / 'geometry'
    geometry_dir.mkdir(exist_ok=True)
    geometry_files = sorted(Path(static_dir).glob("*STATIC_*.h5"))
    crs = get_raster_crs(stitched_ifg_paths[0])
    epsg = crs.to_epsg()
    out_bounds = get_raster_bounds(stitched_ifg_paths[0])
    if water_mask_file is not None:
        if not Path(inps.water_mask_file).exists():
            water_mask_file = create_external_files(water_mask_file,
               stitched_ifg_paths[0], out_bounds, crs, geometry_dir,
               maskfile=True)
    #
    # generate DEM file
    if dem_file is not None:
        if not Path(inps.dem_file).exists():
            dem_file = create_external_files(dem_file,
               stitched_ifg_paths[0], out_bounds, crs, geometry_dir,
               demfile=True)
    #
    # generate geometry files
    frame_geometry_files = prepare_geometry(
        geometry_dir=geometry_dir,
        geo_files=geometry_files,
        matching_file=stitched_ifg_paths[0],
        dem_file=dem_file,
        epsg=epsg,
        out_bounds=out_bounds,
        strides=cfg_obj.output_options.strides
    )
    #
    # unwrap stitched ifg
    row_looks, col_looks = cfg_obj.phase_linking.half_window.to_looks()
    nlooks = row_looks * col_looks
    # set unw options
    unwrap_options = cfg_obj.unwrap_options
    unwrap_options._directory = stitched_ifg_path
    unwrap_options.ntiles = ntiles
    unwrap_options.n_parallel_jobs = n_parallel_jobs
    unwrap_options.n_parallel_tiles = n_parallel_tiles
    unwrap_options.unwrap_method = unwrap_method
    unwrap_options.zero_where_masked = True
    unwrapped_paths, conncomp_paths = unwrapping.run(
        ifg_file_list=stitched_ifg_paths,
        cor_file_list=stitched_cor_paths,
        nlooks=nlooks,
        unwrap_options=unwrap_options,
        temporal_coherence_file=stitched_temp_coh_file,
        mask_file=water_mask_file,
        add_overviews=True
    )
    #
    # go through final time-series stage
    inverted_phase_paths = timeseries.run(
        unwrapped_paths=unwrapped_paths,
        conncomp_paths=conncomp_paths,
        corr_paths=stitched_cor_paths,
        condition_file=stitched_amp_dispersion_file,
        condition=CallFunc.MAX,
        output_dir=stitched_ifg_path / 'timeseries',
        run_velocity=True,
        correlation_threshold=correlation_threshold,
        num_threads=threads_per_worker
    )

    return


if __name__ == '__main__':
    inp = cmd_line_parse()
    access_cslcs(inp)
