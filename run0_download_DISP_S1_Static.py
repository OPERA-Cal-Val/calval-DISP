#!/usr/bin/env python3
import argparse
import os, json
from datetime import datetime
from pathlib import Path

from datetime import datetime as dt

import asf_search as asf
from opera_utils.geometry import stitch_geometry_layers
from opera_utils.download import L2Product

from concurrent.futures import ThreadPoolExecutor, as_completed

import boto3
import botocore

import requests
from io import BytesIO
import zipfile

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

def createParser(iargs = None):
    '''Commandline input parser'''
    parser = argparse.ArgumentParser(description='Downloading OPERA DISP-S1 from AWS S3 bucket and static layer files from ASF')
    parser.add_argument("--frameID", 
                        required=True, type=str, help='frameID of DISP-S1 to download (e.g., 33039)')
    parser.add_argument("--dispDir",
                        default='outputs', type=str, help='directory to download DISP-S1 (default: outputs)')
    parser.add_argument("--startDate", 
                        default='20160101', type=str, help='start date of DISP-S1 (default: 20160101)')
    parser.add_argument("--endDate", 
                        default=dt.today().strftime('%Y%m%d'), type=str, help='end date of DISP-S1 (default: today)')
    parser.add_argument("--nWorkers",
                        default=5, type=int, help='number of simultaenous downloads from AWS S3 bucket (default: 5)')
    parser.add_argument("--staticDir",
                        default='static_lyrs', type=str, help='directory to store static layer files (default: static_lyrs)')
    parser.add_argument("--geomDir",
                        default='geometry', type=str, help='directory to store geometry files from static layers (default: geometry)')
    return parser.parse_args(args=iargs)

def download_file(bucket_name, file_key, local_path):
    ''' download files from S3 bucket '''
    s3 = boto3.client('s3', config=botocore.client.Config(signature_version=botocore.UNSIGNED))
    if not os.path.exists(local_path):
        s3.download_file(bucket_name, file_key, local_path)
        print(f"File downloaded to {local_path}")
    else:
        print(f'{local_path} already exists')
    
def list_s3_directories(bucket_name, directory_name, keyword=None):
    ''' listing directories in bucket '''
    s3 = boto3.client('s3', config=botocore.client.Config(signature_version=botocore.UNSIGNED))
    paginator = s3.get_paginator('list_objects_v2')
    prefix = directory_name if directory_name.endswith('/') else directory_name + '/'
    
    directories = set()
    
    for page in paginator.paginate(Bucket=bucket_name, Prefix=prefix, Delimiter='/'):
        for prefix in page.get('CommonPrefixes', []):
            dir_name = prefix['Prefix']
            if keyword is None or keyword.lower() in dir_name.lower():
                directories.add(dir_name)
    
    return sorted(directories)

def get_key(s):
    return '_'.join(s.split('_')[:-2])  # finding last two separated by underscores

def parse_date(date_string):
    return datetime.strptime(date_string, "%Y%m%d")

def filter_list_by_date_range(list_, start_date, end_date):
    ''' filtered based on start and end date '''

    start = parse_date(start_date)
    end = parse_date(end_date)
    
    filtered_list = []
    for item in list_:
        item_start = parse_date(item[30:38])
        item_end = parse_date(item[47:55])

        if (start <= item_start <= end) and (start <= item_end <= end):
            filtered_list.append(item)
    return filtered_list

def main(inps):
    frameID = inps.frameID
    frameID = frameID.zfill(5)	# force frameID to have 5 digit number as string 
    dispDir = inps.dispDir
    os.makedirs(dispDir, exist_ok='True')
    startDate = inps.startDate
    endDate = inps.endDate
    nWorkers = inps.nWorkers
    staticDir = inps.staticDir
    os.makedirs(staticDir, exist_ok='True')
    geomDir = inps.geomDir
    os.makedirs(geomDir, exist_ok='True')

    bucket_name = 'opera-pst-rs-pop1'       # aws S3 bucket of PST
    directory_name = f'products/DISP_S1'    # directory name where DISP-S1s locate

    print('S3 bucket name: ', bucket_name)
    print('DISP_S1 directory name in bucket: ', directory_name)

    keyword = 'F' + frameID
    subdirectories = list_s3_directories(bucket_name, directory_name, keyword)  # search by frame ID
    list_disp = [ dir.split('/')[-2] for dir in subdirectories]
    list_disp = sorted(list_disp)

    unique_dict = {get_key(x): x for x in list_disp}
    list_disp = list(unique_dict.values())

    list_disp = filter_list_by_date_range(list_disp, startDate, endDate)       # filter by dates

    print('number of DISP-S1 to download: ', len(list_disp))

    # Concurrent downloading of DISP-S1 nc files
    with ThreadPoolExecutor(max_workers=nWorkers) as executor:
        future_to_file = {executor.submit(download_file, bucket_name, f'products/DISP_S1/{select_disp}/{select_disp}.nc', f'{dispDir}/{select_disp}.nc'): select_disp for select_disp in list_disp}
        for future in as_completed(future_to_file):
            select_disp = future_to_file[future]
            try:
                future.result()
            except Exception as exc:
                print(f'{select_disp} generated an exception: {exc}')

    # burst_ids = opera_utils.get_burst_ids_for_frame(int(frameID))     # forced to download zip file to cache

    ## Access json matching bursts to frame IDs without downloading
    # URL of the ZIP file containing the JSON file
    repo_zip_url = 'https://github.com/opera-adt/burst_db/releases/download/v0.5.0/opera-s1-disp-0.5.0-frame-to-burst.json.zip'

    # Access the ZIP file
    response = requests.get(repo_zip_url)
    zip_data = BytesIO(response.content)

    # Extract the JSON file from the ZIP archive
    with zipfile.ZipFile(zip_data, 'r') as zip_ref:
        # Assuming your JSON file is named 'data.json' within the ZIP
        json_data = zip_ref.read('opera-s1-disp-0.5.0-frame-to-burst.json') 

    # Load the JSON data
    data = json.loads(json_data.decode('utf-8')) # ['features']
    burst_ids = data['data'][frameID.lstrip('0')]['burst_ids']  # list of burst IDs within one frame ID

    # search CLSC Static Layer files
    product = L2Product.CSLC_STATIC

    results = asf.search(
        operaBurstID=list(burst_ids),
        processingLevel=product.value,
    )

    results.download(path=staticDir, processes=5)    # downloading static layers with simultaneous downloads

    list_static_files = [ Path(f'{staticDir}/{results[ii].properties["fileName"]}') for ii in range(len(results)) ] 

    print('number of static layer files to download: ', len(results))
    print(list_static_files)

    # generating los_east.tif and los_north.tif from downloaded static layers
    output_files = stitch_geometry_layers(list_static_files, output_dir=geomDir)

    print('Done')

if __name__ == '__main__':
    # load arguments from command line
    inps = createParser()

    print("==================================================================")
    print("        Downloading DISP-S1 and static layer files")
    print("==================================================================")
    
    # Run the main function
    main(inps)
