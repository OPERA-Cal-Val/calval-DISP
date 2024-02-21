#!/usr/bin/env python3
"""Tools for accessing PST CSLC bursts and executing dolphin workflow"""

import argparse
import ast
import datetime
import subprocess
import time
from pathlib import Path

import boto3
from botocore import UNSIGNED
from botocore.config import Config

import geopandas as gpd

import pandas as pd

from shapely.geometry import box

from tqdm import tqdm

from dolphin.io import get_raster_bounds, get_raster_crs

from dolphin.utils import prepare_geometry

from dolphin.workflows import _cli_config as dconfig
from dolphin.workflows import _cli_run as drun
from dolphin.workflows import stitching_bursts, unwrapping, config


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

    return parser


def cmd_line_parse(iargs=None):
    """
        Parser wrapper
    """
    parser = create_parser()
    user_inp = parser.parse_args(args=iargs)

    return user_inp


def download_whole_file(url, local_filename, verbose=False):
    """
        Download specified PST CSLCs
    """
    print("Bulk download of whole file:")
    t0 = time.time()
    # set s3 inputs
    s3_bucket_name = url.split("/")[2]
    s3_path = url.split(s3_bucket_name + '/')[-1]
    # only proceed if corresponding local file does not exist
    if not local_filename.exists():
        s3 = boto3.client('s3', config=Config(signature_version=UNSIGNED))
        # download file
        s3.download_file(s3_bucket_name, s3_path, local_filename)
        if verbose:
            print(f"Took {time.time() - t0:.3f} seconds for bulk download")
    return


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


def access_cslcs(inps=None):
    """
        Custom Dolphin workflow to access and leverage PST CSLCs
    """
    # Initialize variables and workspace
    csv_path = Path(__file__).parent / 'data'
    csv_path = csv_path / \
               'getallbursts_table_validation_bursts_target_v1.0.csv'
    start_date = inps.start_date
    end_date = inps.end_date
    fixed_or_range = inps.fixed_or_range
    frameid_number = inps.frame_id
    orbit_pass = inps.orbit_pass
    area_of_interest = inps.area_of_interest
    strides = inps.strides.split()
    threads_per_worker = inps.threads_per_worker
    n_parallel_jobs = inps.n_parallel_jobs
    n_parallel_tiles = inps.n_parallel_tiles
    ntiles = (inps.ntiles, inps.ntiles)
    if orbit_pass:
        if orbit_pass[:3].lower() == 'asc':
            orbit_pass = 'Ascending'
        if orbit_pass[:3].lower() == 'des':
            orbit_pass = 'Descending'

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

    # convert the string to a datetime object
    date_format = "%Y%m%d"
    if start_date and end_date:
        start_datetime = datetime.datetime.strptime(start_date, date_format)
        end_datetime = datetime.datetime.strptime(end_date, date_format)

    #
    # Open and setup dataframe
    df = pd.read_csv(csv_path)
    # parse geometry
    df = gpd.GeoDataFrame(
        df.loc[:, [c for c in df.columns if c != 'geometry']],
        geometry=gpd.GeoSeries.from_wkt(df['geometry'])
        )
    # parse datetime
    df['date_datetime'] = pd.to_datetime(df['date'], format=date_format)
    # read frame IDs as list
    df['frame_ids'] = df['frame_ids'].apply(lambda x: ast.literal_eval(x))

    #
    # Apply filters to dataframe
    # filter by geometry
    if area_of_interest:
        df = df[df['geometry'].apply(lambda x:
                area_of_interest.intersects(x))]

    # filter by orbit direction
    if orbit_pass:
        df = df[df['orbit_pass_direction'].apply(lambda x: orbit_pass == x)]

    # If specified, only capture input dates which match input
    # as opposed to treating them as a range
    if start_date and end_date:
        if fixed_or_range:
            # Filter dataframe by specified date range
            df = df[(df['date_datetime'] == start_datetime)
                    | (df['date_datetime'] == end_datetime)]
        else:
            df = df[(df['date_datetime'] >= start_datetime)
                    & (df['date_datetime'] <= end_datetime)]

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
    if frameid_number:
        df = df[df['frame_ids'].apply(lambda x: frameid_number in x)]

    # Reset the index of the filtered DataFrame
    df = df.reset_index(drop=True)

    #
    # Produce burstwise IFGs over filtered subset
    # run dolphin by burst
    list_of_bursts = list(set(df['burst_id'].to_list()))
    for burst in tqdm(list_of_bursts, desc="Processing burst:"):
        print(f"{burst}")
        # set output path for burst
        dolphin_dir_burst = dolphin_dir.joinpath(f'{burst}')
        if not dolphin_dir_burst.exists():
            dolphin_dir_burst.mkdir(parents=True, exist_ok=True)
        # filter by burst
        filtered_df = df[df['burst_id'] == burst]
        filtered_df = filtered_df.reset_index(drop=True)
        cslc_txt = out_dir.joinpath(f'CSLC_list_{burst}.txt')
        static_txt = out_dir.joinpath(f'static_list_{burst}.txt')
        for index, row in tqdm(filtered_df.iterrows(),
                               desc="Downloading scene"):
            print(f"{row['date']}")
            # download CSLC
            url = row['cslc_url']
            local_filename = url.split("/")[-1]
            local_filename = cslc_dir.joinpath(local_filename)
            download_whole_file(url, local_filename)
            # record downloads in respective text-files
            target_prefix = 'OPERA_L2_CSLC-S1_T'
            create_inpt_txtfile(cslc_dir, target_prefix, cslc_txt)
            if index == 0:
                # download static layer just for the first date
                url = row['cslc_static_url']
                local_filename = url.split("/")[-1]
                local_filename = static_dir.joinpath(local_filename)
                download_whole_file(url, local_filename)
                # record downloads in respective text-files
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
                              mask_file=None
                              no_unwrap=True)
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

    #
    # Stitch all ifgs by burst
    # set static paths for dolphin
    ifg_paths = list(dolphin_dir.glob('t*/interferograms/*.vrt'))
    coh_paths = list(dolphin_dir.glob('t*/linked_phase/' +
                     'temporal_coherence_average_*.tif'))
    # only access looked file if applicable
    if strides == ['1', '1']:
        ps_file_list = list(dolphin_dir.glob('t*/PS/ps_pixels.tif'))
    else:
        ps_file_list = list(dolphin_dir.glob('t*/PS/ps_pixels_looked.tif'))
    cfg_obj = config.DisplacementWorkflow.from_yaml(yml_file)
    (
        stitched_ifg_paths,
        stitched_cor_paths,
        _,
        _,
    ) = stitching_bursts.run(
        ifg_file_list=ifg_paths,
        temp_coh_file_list=coh_paths,
        ps_file_list=ps_file_list,
        stitched_ifg_dir=stitched_ifg_path,
        output_options=cfg_obj.output_options,
        file_date_fmt='%Y%m%d',
        corr_window_size=(11, 11),
    )
    #
    # generate geometry files
    dem_file = None #!#
    geometry_dir = stitched_ifg_path / 'geometry'
    geometry_dir.mkdir(exist_ok=True)
    geometry_files = sorted(Path(static_dir).glob("*STATIC_*.h5"))
    crs = get_raster_crs(stitched_ifg_paths[0])
    epsg = crs.to_epsg()
    out_bounds = get_raster_bounds(stitched_ifg_paths[0])
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
    # unwrap stitched unw
    row_looks, col_looks = cfg_obj.phase_linking.half_window.to_looks()
    nlooks = row_looks * col_looks
    # set unw options
    unwrap_options = cfg_obj.unwrap_options
    unwrap_options._directory = stitched_ifg_path
    unwrap_options.ntiles = ntiles
    unwrap_options.n_parallel_jobs = n_parallel_jobs
    unwrap_options.n_parallel_tiles = n_parallel_tiles
    unwrapping.run(
        ifg_file_list=stitched_ifg_paths,
        cor_file_list=stitched_cor_paths,
        nlooks=nlooks,
        unwrap_options=unwrap_options,
        mask_file=cfg_obj.mask_file,
    )

    return


if __name__ == '__main__':
    inp = cmd_line_parse()
    access_cslcs(inp)
