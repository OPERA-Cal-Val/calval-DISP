#!/usr/bin/env python
'''
Stitch DISP products and generate basic report for qualitative checks
'''

import argparse
import glob
import os
import re
import sys
from collections import defaultdict
from os import fspath
from pathlib import Path

import geopandas as gpd
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import osgeo
from pyproj import CRS, Proj, Transformer

from matplotlib.backends.backend_pdf import PdfPages

# Add the src directory to sys.path
sys.path.append(str(Path(__file__).parent / 'src'))

from pst_dolphin_utils import get_raster_crs
from seq_stitch_disp import product_stitch_sequential

def create_parser():
    """
        Run dolphin workflow with specified input parameters
    """
    parser = argparse.ArgumentParser(description='Sort and stitch '
                                     'DISP products')
    parser.add_argument(
        '-f',
        '--file-glob',
        type=str,
        default='./products/*.nc',
        help=f'path pattern of DISP products.'
    )
    parser.add_argument(
        '-ao', '--areaofinterest',
        dest='area_of_interest', type=str,
        default=None,
        help='Specify area of interest as: '
            '1. path to a valid shapefile, or '
            '2. S N W E coordinates, ie ISCE convention'
    )
    parser.add_argument(
        '-o', '--outdir',
        dest='out_dir', type=str,
        default='./',
        help='Specify output directory')
    parser.add_argument(
        '-of', '--outputformat',
        dest='output_format', type=str,
        default='ENVI',
        help='Specify output raster format')

    return parser


def cmd_line_parse(iargs=None):
    """
        Parser wrapper
    """
    parser = create_parser()
    user_inp = parser.parse_args(args=iargs)

    return user_inp


def disp_dict(filenames):
    # Dictionary to store lists of filenames by date pair
    date_pairs_dict = defaultdict(list)

    # Regular expression patterns to extract the date pair and the F number
    date_pattern = re.compile(r"_(\d{8}T\d{6}Z)_(\d{8}T\d{6}Z)")
    f_number_pattern = re.compile(r"_F(\d{5})_")

    # Temporary storage for checking continuity
    temp_storage = defaultdict(list)

    # Extract the date pairs and F numbers
    for filename in filenames:
        date_match = date_pattern.search(filename)
        f_number_match = f_number_pattern.search(filename)
    
        if date_match and f_number_match:
            date_pair = date_match.group(1)[:8] + '_' + date_match.group(2)[:8]
            f_number = int(f_number_match.group(1))
            temp_storage[date_pair].append((f_number, filename))

    # Check for continuity and add to the dictionary if continuous
    for date_pair, values in temp_storage.items():
        f_numbers = sorted([f_number for f_number, _ in values])
    
        # Check if the f_numbers are continuous
        if all(f_numbers[i+1] - f_numbers[i] == 1 for i in range(len(f_numbers) - 1)):
            for _, filename in values:
                date_pairs_dict[date_pair].append(filename)

    # Sort the dictionary by the date pair keys
    sorted_date_pairs = dict(sorted(date_pairs_dict.items()))

    # Print the sorted dictionary
    print(sorted_date_pairs)

    return sorted_date_pairs


def main(iargs=None):
    inps = cmd_line_parse(iargs)

    file_glob = inps.file_glob.split(',')
    area_of_interest = inps.area_of_interest
    output_format = inps.output_format
    out_dir = inps.out_dir
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    if area_of_interest is not None:
        if Path(area_of_interest).is_file():
            area_of_interest = gpd.read_file(area_of_interest)
            area_of_interest = area_of_interest['geometry'][0].bounds
            # convert to S,N,W,E
            area_of_interest = [area_of_interest[1],
                                area_of_interest[3],
                                area_of_interest[0],
                                area_of_interest[2]]
        else:
            area_of_interest = area_of_interest.split()
            area_of_interest = [float(i) for i in area_of_interest]

    unw_pdf = os.path.join(out_dir, 'unw_pngs.pdf')
    if os.path.exists(unw_pdf):
        raise Exception(f'Final output file {unw_pdf} already exists in '
                        f'specified output directory {out_dir}. Please ' 
                        'specify a new output directory or delete the PDF.')

    unw_dir = os.path.join(out_dir, 'unwrappedPhase')
    if not os.path.exists(unw_dir):
        os.mkdir(unw_dir)

    conncomp_dir = os.path.join(out_dir, 'connectedComponents')
    if not os.path.exists(conncomp_dir):
        os.mkdir(conncomp_dir)

    input_files = []
    for i in file_glob:
        input_files.extend(glob.glob(i))
    input_files = sorted(input_files)
    input_files = disp_dict(input_files)

    # loop through each pair
    epsg_code = []
    unw_pngs = []
    for pair_name in input_files:
        # get iterative list of products
        pair_list = input_files[pair_name]
        phs_files = \
            [f'NETCDF:"{i}":unwrapped_phase' for i in pair_list]
        conn_files = \
            [f'NETCDF:"{i}":connected_component_labels' for i in pair_list]

        # get median epsg value for first set of pairs
        if epsg_code == []:
            for i in phs_files:
                epsg_code.append(get_raster_crs(i).to_epsg())
            epsg_code = int(np.median(epsg_code))
            # capture bounds once, if defined
            if area_of_interest is not None:
                # convert bounds to native epsg
                south, north, west, east = area_of_interest
                transformer = Transformer.from_crs('EPSG:4326',
                    f'EPSG:{epsg_code}', always_xy=True)
                transformed_bbox = [transformer.transform(lon, lat) for lon, lat in
                    [(west, south), (east, south), (west, north), (east, north)]]
                west_convtd = min(pt[0] for pt in transformed_bbox)
                south_convtd = min(pt[1] for pt in transformed_bbox)
                east_convtd = max(pt[0] for pt in transformed_bbox)
                north_convtd = max(pt[1] for pt in transformed_bbox)
                bounds = [west_convtd, south_convtd, east_convtd, north_convtd]
            else:
                bounds = None

        # automatically set other variables
        outFilePhs = os.path.join(unw_dir, pair_name)
        outFileConnComp = os.path.join(conncomp_dir, pair_name)

        # stitching
        unw_pngs.append(outFilePhs + '.png')
        product_stitch_sequential(
            phs_files, conn_files,
            arrres=[30., 30.], epsg=epsg_code,
            output_unw=outFilePhs,
            output_conn=outFileConnComp,
            output_format=output_format,
            bounds=bounds,
            range_correction=True, save_fig=True,
            overwrite=True)

    # Create PDF with all images concatenated
    with PdfPages(unw_pdf) as pdf:
        for unw in unw_pngs:
            # Read the image
            img = mpimg.imread(unw)
            # Create a figure and disable the axis
            fig, ax = plt.subplots()
            ax.axis('off')
            # Display the image
            ax.imshow(img)
            # Save the current figure to the PDF
            pdf.savefig(fig)
            # Close the figure
            plt.close(fig)

if __name__ == "__main__":
    main(sys.argv[1:])
