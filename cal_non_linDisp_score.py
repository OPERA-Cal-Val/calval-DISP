#!/usr/bin/env python3
import argparse
import os
import numpy as np
import time
import pandas as pd
import matplotlib.pyplot as plt
import psutil

import h5py

from typing import Tuple
from mintpy.utils import writefile
import concurrent.futures
from tqdm import tqdm

import warnings
warnings.filterwarnings("ignore")

ncpus = len(psutil.Process().cpu_affinity())    # number of available CPUs

# Define the function to be executed in parallel
def calculate_window_std(start_row: int, end_row: int, start_col: int, end_col: int,depth: int, ifg_stack_file: str, window_size: int) -> Tuple[int, int, int, int, int, np.ndarray]:
    ''' 
    calculating standard deviation in temporally windowed data
    '''

    with h5py.File(ifg_stack_file, 'r') as f:
        window_data = f['unwrapPhase'][depth:(depth+window_size), start_row:end_row, start_col:end_col] * (-wavelength/(4*np.pi))      # converting to meter
        window_data = np.transpose(window_data, (1,2,0))      # reshaping an array to have number of pairs to last axis
    
    window_data = window_data.reshape(-1,window_size)   # converting from 3D to 2D
    rolling_mean = pd.DataFrame(window_data).rolling(window=window_size, min_periods=1, axis=1).mean().values   # rolling mean in time
    window_std = np.std(rolling_mean, axis=-1)  # 1D
    window_std = window_std.reshape(end_row-start_row, end_col-start_col)   # converting to 2D

    return start_row, end_row, start_col, end_col, depth, window_std

def createParser(iargs = None):
    '''Commandline input parser'''
    parser = argparse.ArgumentParser(description='Generating non-linear displacement score map')
    parser.add_argument("--mintpyDir",
                        default='mintpy_output', type=str, help='directory containing MintPy hdf files (default: mintpy_output)')
    parser.add_argument("--ifgStack", 
                        default='inputs/ifgramStack.h5', type=str, help='location ifg stack hdf file inside mintpyDir (default: inputs/ifgramStack.h5)')
    parser.add_argument("--scoreMap", 
                        default='nonDispScore.h5', type=str, help='output: normalized non-linear displacement temporal score map (default: nonDispScore.h5)')
    parser.add_argument("--winSize", 
                        default=11, type=int, help='temporal window size for calculating score (default: 11)')    
    parser.add_argument("--nProcs", 
                        default=ncpus, type=int, help='number of used cpu-cores (default: number of all available CPUs)')   
    parser.add_argument("--nPatch", 
                        default=20, type=int, help='number of patches in each row and col for multi processing (default: 20)')
    parser.add_argument("--nDepth", help='number of pairs for calculating non-linear displacement score map (defaut: None, meaning all data are used)')
    return parser.parse_args(args=iargs)

def main(inps):
    mintpyDir = inps.mintpyDir
    ifgStack = inps.ifgStack
    scoreMap = inps.scoreMap
    winSize = inps.winSize
    nProcs = inps.nProcs
    nPatch = inps.nPatch
    nDepth = inps.nDepth

    ifgStack = f'{mintpyDir}/{ifgStack}'
    scoreMap = f'{mintpyDir}/{scoreMap}'

    print(f'number of all available CPU cores: {ncpus}')
    print(f'number of CPU cores to be used: {nProcs} \n')

    assert nProcs <= ncpus  # input no. of cpus should be smaller than no. of available cpus

    if os.path.exists(scoreMap):
        print(f'{scoreMap} already exists')
    else:
        with h5py.File(ifgStack, 'r') as f:
            atr = dict(f.attrs)
            global wavelength
            wavelength = float(atr['WAVELENGTH'])
            num_ifgs , row, col = f['unwrapPhase'].shape
        
        if nDepth:  # if number of depth is set
            num_ifgs = int(nDepth)      # replacing number of interferograms with input number

        step_col = col // nPatch
        step_row = row // nPatch

        list_col = list(range(0, col, step_col))
        list_row = list(range(0, row, step_row))

        ncol = len(list_col)
        nrow = len(list_row)

        # constructing parameters in patches for multi-processing
        params = []

        for ind_col, start_col in enumerate(list_col):
            for ind_row, start_row in enumerate(list_row):

                if ind_col == (ncol - 1):
                    end_col = col
                else:
                    end_col = list_col[ind_col + 1]

                if ind_row == (nrow - 1):
                    end_row = row
                else:
                    end_row = list_row[ind_row + 1]

                for depth in range(num_ifgs - winSize+ 1):
                    params.append((start_row, end_row, start_col, end_col, depth, ifgStack, winSize))

        print('number of multi-processors: ', len(params))

        variability_indices = np.zeros((row, col, num_ifgs - winSize + 1), dtype=np.float32)    # allocating space 

        st_time = time.time()

        # excuting multi-processing with patches
        with concurrent.futures.ProcessPoolExecutor(max_workers=nProcs) as executor:
            futures = [executor.submit(calculate_window_std, *param) for param in params]
            
            for future in tqdm(concurrent.futures.as_completed(futures), total=len(params), desc="Processing"):
                result = future.result()
                start_row, end_row, start_col, end_col, depth, window_std = result
                variability_indices[start_row:end_row, start_col:end_col ,depth] = window_std

        end_time = time.time()
        time_taken = np.round((end_time - st_time)/60.,2)
        print(f'{time_taken} min taken for multi-core non-linear masking processing')

        sigma = np.std(variability_indices,axis=-1)
        variability_scores = np.sum(variability_indices < sigma[..., np.newaxis], axis=2) / (num_ifgs - winSize + 1)     # ratio: normalized between 0 and 1

        fig, ax = plt.subplots(figsize=(12, 12))
        im = ax.imshow(variability_scores, cmap='viridis', interpolation='none')
        fig.colorbar(im, ax=ax, shrink=0.2)
        ax.set_title('Temporal variability score map')
        fig.savefig(f'{mintpyDir}/score_map.png', dpi=300, bbox_inches='tight', pad_inches = 0)

        var_atr = {}
        var_atr['WIDTH'] = col
        var_atr['LENGTH'] = row
        var_atr['FILE_TYPE'] = 'mask'
        var_atr['DATA_TYPE'] = np.float32
        writefile.write(variability_scores, out_file = scoreMap, metadata=var_atr)  # writing variability score file

    print('Done')

if __name__ == '__main__':
    # load arguments from command line
    inps = createParser()

    print("==================================================================")
    print("        Generating the non-linear displacement score map")
    print("==================================================================")
    
    # Run the main function
    main(inps)
