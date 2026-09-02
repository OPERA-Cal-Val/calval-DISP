#!/usr/bin/env python
"""
run2b_stitch_track_frames.py

Stitch multiple spatially continuous InSAR velocity frames across a single track.

This utility is intended to be run after generating velocity and timeseries 
outputs (via run2_prep_mintpy_opera.py). 

USER RESPONSIBILITY:
- All provided frames MUST belong to the same track.
- Frames MUST be spatially continuous. 
- The script assumes a strict directory structure where parent directories are 
  named with their OPERA Frame ID (e.g., F11115, F11116)[cite: 2], and consecutive ID 
  numbers denote spatially overlapping frames.
"""

import argparse
import glob
import os
import re
import numpy as np
import rasterio
from rasterio.merge import merge
from rasterio.windows import from_bounds, Window
from osgeo import gdal
from mintpy.cli import save_gdal


def parse_args():
    parser = argparse.ArgumentParser(
        description="Stitch multiple spatially continuous InSAR velocity frames across a single track.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        '-i', '--input-files', type=str, required=True, nargs='+',
        help="Input velocity.tif files (wildcards supported). \n"
             "Directories must follow the FXXXXX naming convention[cite: 2]."
    )
    parser.add_argument(
        '-o', '--output-dir', type=str, required=True,
        help="Destination directory for the final stitched products."
    )
    parser.add_argument(
        '-r', '--ref-lat-lon', type=str, required=True,
        help="Final spatial reference point passed as 'lat lon' (e.g., '40.32 -123.13')[cite: 2]."
    )
    parser.add_argument(
        '--skip-mask', action='store_true',
        help="Skip applying the 90% reliability threshold mask to the velocity files."
    )
    return parser.parse_args()


def extract_and_reproject(track_path, out_dir, frame_id, skip_mask=False):
    """
    Extracts MintPy HDF5 datasets to TIF, applies an optional reliability mask, 
    and reprojects the files to WGS84 on a common snapped grid.
    """
    # 1. Check and extract H5 files
    h5_files = {
        'temporalCoherence': f'{track_path}/avg_lyrs/temporalCoherence.h5',
        'mask': f'{track_path}/recommended_mask_90thresh.h5'
    }
    
    tif_paths = {}
    for name, h5_path in h5_files.items():
        if not os.path.exists(h5_path):
            h5_path = h5_path.replace('avg_lyrs/', '')  # fallback
            
        tif_path = f'{track_path}/{name}.tif'
        if not os.path.exists(tif_path):
            print(f"Extracting {name} to TIF for {frame_id}...")
            iargs = [h5_path, '-d', name, '-o', tif_path]
            
            # ONLY apply zero-masking to float datasets like coherence.
            if name != 'mask':
                iargs.append('--zm')
                
            save_gdal.main(iargs)
        tif_paths[name] = tif_path

    # 2. Identify velocity files and apply mask if requested
    vel_files = {
        'vel': f'{track_path}/velocity.tif',
        'vel_corr': f'{track_path}/velocity_demErr_iono_SET.tif'
    }
    
    files_to_warp = {'coh': tif_paths['temporalCoherence']}
    
    if skip_mask:
        print(f"Skipping mask application for {frame_id}...")
        for key, v_path in vel_files.items():
            if os.path.exists(v_path):
                files_to_warp[key] = v_path
    else:
        mask_data = None
        for key, v_path in vel_files.items():
            if not os.path.exists(v_path):
                continue
            out_mask_path = v_path.replace('.tif', '_msk.tif')
            
            if not os.path.exists(out_mask_path):
                print(f"Masking {key} for {frame_id}...")
                if mask_data is None:
                    with rasterio.open(tif_paths['mask']) as src_mask:
                        mask_data = src_mask.read(1)
                        
                with rasterio.open(v_path) as src_v:
                    meta = src_v.meta.copy()
                    v_data = src_v.read(1)
                    
                    # Apply mask and enforce nodata tag
                    v_data[mask_data == 0] = np.nan
                    meta.update({"nodata": np.nan})
                    
                    with rasterio.open(out_mask_path, 'w', **meta) as dst_v:
                        dst_v.write(v_data, 1)
                        
            files_to_warp[key] = out_mask_path

    # 3. Reproject everything to WGS84 on a common snapped grid[cite: 2]
    warp_opts = gdal.WarpOptions(
        dstSRS='EPSG:4326', format='GTiff',
        xRes=0.00027, yRes=0.00027, targetAlignedPixels=True,
        creationOptions=['COMPRESS=LZW'],
        dstNodata=np.nan  # Force empty borders to be NaN
    )
    
    wgs84_files = {}
    for key, fpath in files_to_warp.items():
        out_wgs84 = os.path.join(out_dir, f"{frame_id}_{key}_wgs84.tif")
        if not os.path.exists(out_wgs84):
            print(f"Reprojecting {fpath} to WGS84...")
            gdal.Warp(out_wgs84, fpath, options=warp_opts)
        wgs84_files[key] = out_wgs84
        
    return wgs84_files


def get_tie_point_and_offset(f1_vel, f2_vel, f1_coh, f2_coh):
    """
    Finds the overlapping region between two consecutive frames, 
    evaluates coherence, and calculates the robust median offset.
    """
    with rasterio.open(f1_vel) as src1, rasterio.open(f2_vel) as src2, \
         rasterio.open(f1_coh) as coh1, rasterio.open(f2_coh) as coh2:
             
        # Find intersection of ALL FOUR datasets to be absolutely safe
        left = max(src1.bounds.left, src2.bounds.left, coh1.bounds.left, coh2.bounds.left)
        bottom = max(src1.bounds.bottom, src2.bounds.bottom, coh1.bounds.bottom, coh2.bounds.bottom)
        right = min(src1.bounds.right, src2.bounds.right, coh1.bounds.right, coh2.bounds.right)
        top = min(src1.bounds.top, src2.bounds.top, coh1.bounds.top, coh2.bounds.top)
        
        if left >= right or bottom >= top:
            print("\nWARNING: Bounding boxes do not overlap at all!")
            return None, None, 0.0

        # 1. Generate the base window from src2 (the new frame)
        base_win = from_bounds(left, bottom, right, top, src2.transform).round_offsets().round_shape()
        
        # Read src2 array to get the EXACT pixel dimensions of the overlapping region
        v2_arr = src2.read(1, window=base_win, boundless=True, fill_value=np.nan)
        h, w = v2_arr.shape

        # 2. Force all other datasets to read the EXACT same pixel dimensions 
        def get_matched_window(dataset):
            win = from_bounds(left, bottom, right, top, dataset.transform).round_offsets()
            return Window(win.col_off, win.row_off, w, h)

        win1_v = get_matched_window(src1)
        win1_c = get_matched_window(coh1)
        win2_c = get_matched_window(coh2)
        
        v1_arr = src1.read(1, window=win1_v, boundless=True, fill_value=np.nan)
        c1_arr = coh1.read(1, window=win1_c, boundless=True, fill_value=np.nan)
        c2_arr = coh2.read(1, window=win2_c, boundless=True, fill_value=np.nan)
        
        # Find valid pixels in both
        valid = (
            np.isfinite(v1_arr) & (v1_arr != 0) &
            np.isfinite(v2_arr) & (v2_arr != 0) &
            np.isfinite(c1_arr) & np.isfinite(c2_arr)
        )
        
        if not np.any(valid):
            print("\nWARNING: Bounding boxes overlap, but valid data does not!")
            print("Assuming an offset of 0.0 and continuing the stitch...")
            return None, None, 0.0
            
        combined_coh = c1_arr * c2_arr
        combined_coh[~valid] = -1
        
        # Grab the top 20% most coherent pixels across the entire overlap[cite: 2]
        coh_threshold = np.percentile(combined_coh[valid], 80)
        best_pixels = valid & (combined_coh >= coh_threshold)
        
        if not np.any(best_pixels):
            best_pixels = valid
            
        # Calculate the robust median offset using ALL the best pixels
        v1_best = v1_arr[best_pixels]
        v2_best = v2_arr[best_pixels]
        offset = float(np.median(v1_best - v2_best))
        
        # Log the absolute max coherence point[cite: 2]
        max_idx = np.unravel_index(np.argmax(combined_coh), combined_coh.shape)
        x_off = win1_v.col_off + max_idx[1]
        y_off = win1_v.row_off + max_idx[0]
        tie_lon, tie_lat = src1.xy(y_off, x_off)
        
        return tie_lat, tie_lon, offset


def apply_offset_and_stitch(f1_path, f2_path, offset, out_path):
    """Applies the calculated offset to the incoming frame and merges it."""
    with rasterio.open(f2_path) as src2:
        meta = src2.meta.copy()
        v2_arr = src2.read(1)
        valid = np.isfinite(v2_arr) & (v2_arr != 0)
        v2_arr[valid] += offset
        
        # Force all invalid regions to strict NaN for clean merge[cite: 2]
        v2_arr[~valid] = np.nan
        meta.update({"nodata": np.nan})
        
        shifted_f2 = f2_path.replace('.tif', '_shifted.tif')
        with rasterio.open(shifted_f2, 'w', **meta) as dst2:
            dst2.write(v2_arr, 1)
            
    src_files_to_mosaic = [rasterio.open(f1_path), rasterio.open(shifted_f2)]
    
    # nodata=np.nan forces merge to ignore overlapping empty pixels[cite: 2]
    mosaic, out_trans = merge(src_files_to_mosaic, nodata=np.nan)
    
    out_meta = src_files_to_mosaic[0].meta.copy()
    out_meta.update({
        "height": mosaic.shape[1],
        "width": mosaic.shape[2],
        "transform": out_trans,
        "nodata": np.nan
    })
    
    with rasterio.open(out_path, "w", **out_meta) as dest:
        dest.write(mosaic)
        
    for src in src_files_to_mosaic:
        src.close()
    os.remove(shifted_f2)


def final_reference_shift(merged_path, ref_lat, ref_lon):
    """Shifts the final map based on the designated reference coordinate."""
    with rasterio.open(merged_path) as src:
        arr = src.read(1)
        meta = src.meta.copy()
        
        row, col = src.index(ref_lon, ref_lat)
        ref_val = arr[row, col]
        
        if np.isnan(ref_val) or ref_val == 0:
            print(f"Warning: Reference pt {ref_lat},{ref_lon} is NaN/0!")
            
        print(
            f"Referencing map to {ref_lat}, {ref_lon} "
            f"(subtracting {ref_val:.4f})..."
        )
        
        valid = np.isfinite(arr) & (arr != 0)
        arr[valid] -= ref_val
        
    final_path = merged_path.replace('.tif', '_referenced.tif')
    with rasterio.open(final_path, 'w', **meta) as dst:
        dst.write(arr, 1)


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    
    ref_lat, ref_lon = map(float, args.ref_lat_lon.split())
    
    in_files = []
    for f in args.input_files:
        in_files.extend(glob.glob(f))
        
    frames = []
    for fpath in in_files:
        track_path = os.path.dirname(fpath)
        match = re.search(r'F\d{5}', track_path)[cite: 2]
        if match:
            frames.append((match.group(), track_path))
            
    frames.sort(key=lambda x: x[0])
    print(f"Found contiguous frames: {[f[0] for f in frames]}")
    
    if not frames:
        print("No valid frames found. Exiting.")
        return

    final_stitched_name = "+".join([f[0] for f in frames])
    final_stitched_vel = os.path.join(
        args.output_dir, f"frames{final_stitched_name}_vel.tif"
    )
    final_stitched_corr = os.path.join(
        args.output_dir, f"frames{final_stitched_name}_vel_corr.tif"
    )

    # Bypass extraction/stitching if final mosaic exists
    if os.path.exists(final_stitched_vel):
        print(f"\nFinal stitched map {final_stitched_name} already exists.")
        print("Skipping extraction/stitching. Proceeding directly to referencing...")
        current_stitched_vel = final_stitched_vel
        current_stitched_corr = (
            final_stitched_corr if os.path.exists(final_stitched_corr) else None
        )
    else:
        processed_files = {}

        def get_extracted(f_id, f_path):
            if f_id not in processed_files:
                processed_files[f_id] = extract_and_reproject(
                    f_path, args.output_dir, f_id, args.skip_mask
                )
            return processed_files[f_id]

        f1_id, f1_path = frames[0]
        stitched_name = f1_id
        
        current_stitched_vel = os.path.join(args.output_dir, f"{f1_id}_vel_wgs84.tif")
        current_stitched_corr = os.path.join(args.output_dir, f"{f1_id}_vel_corr_wgs84.tif")
        if not os.path.exists(current_stitched_vel):
            f1_files = get_extracted(f1_id, f1_path)
            current_stitched_vel = f1_files['vel']
            current_stitched_corr = f1_files.get('vel_corr')
        elif not os.path.exists(current_stitched_corr):
            current_stitched_corr = None

        for i in range(1, len(frames)):
            f2_id, f2_path = frames[i]
            next_stitched_name = f"{stitched_name}+{f2_id}"
            next_stitched_vel = os.path.join(
                args.output_dir, f"frames{next_stitched_name}_vel.tif"
            )
            next_stitched_corr = os.path.join(
                args.output_dir, f"frames{next_stitched_name}_vel_corr.tif"
            )
            
            if os.path.exists(next_stitched_vel):
                print(f"\nIntermediate stitch {next_stitched_name} already exists. Skipping...")
                stitched_name = next_stitched_name
                current_stitched_vel = next_stitched_vel
                if os.path.exists(next_stitched_corr):
                    current_stitched_corr = next_stitched_corr
                continue
                
            print(f"\nProcessing {f2_id} to stitch with {stitched_name}...")
            f2_files = get_extracted(f2_id, f2_path)
            
            prev_f_id, prev_f_path = frames[i-1]
            prev_f_files = get_extracted(prev_f_id, prev_f_path)

            print(f"Finding tie point between {stitched_name} and {f2_id}...")
            tie_lat, tie_lon, offset = get_tie_point_and_offset(
                current_stitched_vel, 
                f2_files['vel'],
                prev_f_files['coh'],  
                f2_files['coh']
            )
            
            if tie_lat is None:
                print(f"Offset: {offset:.5f} (DEFAULTED TO 0.0 DUE TO NO VALID OVERLAP)")
            else:
                print(
                    f"Tie Point: Lat {tie_lat:.6f}, Lon {tie_lon:.6f} | "
                    f"Offset: {offset:.5f}"
                )
            
            apply_offset_and_stitch(
                current_stitched_vel, f2_files['vel'], 
                offset, next_stitched_vel
            )
            
            if current_stitched_corr and 'vel_corr' in f2_files:
                apply_offset_and_stitch(
                    current_stitched_corr, f2_files['vel_corr'], 
                    offset, next_stitched_corr
                )
                current_stitched_corr = next_stitched_corr
                
            current_stitched_vel = next_stitched_vel
            stitched_name = next_stitched_name

    print("\nApplying final reference point...")
    final_reference_shift(current_stitched_vel, ref_lat, ref_lon)
    if current_stitched_corr:
        final_reference_shift(current_stitched_corr, ref_lat, ref_lon)
        
    print("Stitching and referencing complete!")


if __name__ == "__main__":
    main()
