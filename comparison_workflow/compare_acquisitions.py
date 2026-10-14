#!/usr/bin/env python

import argparse
import os
import glob
import subprocess
import sys
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime, timedelta
from collections import Counter


def parse_args():
    """Set up command line arguments with defaults."""
    parser = argparse.ArgumentParser(
        description="Compare DAAC and SDS OPERA acquisition dates."
    )
    parser.add_argument(
        "--frame",
        type=str,
        default="F11114",
        help="Frame ID to filter and plot (default: F11114)"
    )
    parser.add_argument(
        "--daac_dir",
        type=str,
        default="DAAC_run/products",
        help="Local path to downloaded DAAC .nc files"
    )
    parser.add_argument(
        "--s3_path",
        type=str,
        default="s3://opera-dev-rs-fwd-gmanipon/products/DISP_S1/",
        help="Base S3 path to search"
    )
    parser.add_argument(
        "--sds_local_dir",
        type=str,
        default=".",
        help="Local path to downloaded SDS .nc files for symlinks"
    )
    parser.add_argument(
        "--sds_daac_combo",
        type=str,
        default="forward+DAAC/products",
        help="Output dir for minimum spanning network symlinks"
    )
    return parser.parse_args()


def fetch_s3_pairs(frame_id, s3_path):
    """Fetches .nc file dates and filenames from S3."""
    cmd = (
        f"aws s3 --profile saml-pub ls {s3_path} --recursive "
        f"| grep '{frame_id}.*\\.nc$'"
    )
    print(f"Running command:\n{cmd}\n")
    
    result = subprocess.run(
        cmd, shell=True, capture_output=True, text=True
    )
    
    s3_pairs = []
    s3_file_map = {}
    lines = result.stdout.strip().split('\n')
    
    for line in lines:
        if not line:
            continue
        filename = line.split()[-1].split('/')[-1]
        parts = filename.split('_')
        
        dates = [
            p for p in parts
            if 'T' in p and p.endswith('Z') and len(p) == 16
        ]
        
        if len(dates) >= 2:
            pair = (dates[0][:8], dates[1][:8])
            s3_pairs.append(pair)
            if pair not in s3_file_map:
                s3_file_map[pair] = []
            s3_file_map[pair].append(filename)
            
    return sorted(s3_pairs), s3_file_map


def fetch_sds_local_pairs(sds_dir, frame_id):
    """Fetches .nc file dates and filenames from a local SDS directory."""
    sds_pairs = []
    sds_file_map = {}
    
    # Check for a 'products' subfolder as expected by the symlink structure
    target_dir = os.path.join(sds_dir, "products")
    if not os.path.isdir(target_dir):
        target_dir = sds_dir
        
    search_pattern = os.path.join(target_dir, f"*{frame_id}*.nc")
    filepaths = glob.glob(search_pattern)
    
    if not filepaths:
        return [], {}, target_dir

    for filepath in filepaths:
        filename = os.path.basename(filepath)
        filename_clean = filename.replace('.nc', '')
        parts = filename_clean.split('_')
        
        dates = [
            p for p in parts
            if 'T' in p and p.endswith('Z') and len(p) == 16
        ]
        
        if len(dates) >= 2:
            pair = (dates[0][:8], dates[1][:8])
            if pair not in sds_file_map:
                sds_file_map[pair] = []
            sds_file_map[pair].append(filename)
            
    for pair in sds_file_map.keys():
        sds_pairs.append(pair)
        
    return sorted(sds_pairs), sds_file_map, target_dir


def parse_daac_directory(daac_dir, frame_id):
    """Reads local .nc files to extract date pairs for a frame."""
    daac_pairs = []
    daac_file_map = {}
    search_pattern = os.path.join(daac_dir, f"*{frame_id}*.nc")
    filepaths = glob.glob(search_pattern)
    
    if not filepaths:
        print(f"Warning: No valid .nc files found in {daac_dir}")
        return daac_pairs, daac_file_map

    for filepath in filepaths:
        filename = os.path.basename(filepath)
        filename_clean = filename.replace('.nc', '')
        parts = filename_clean.split('_')
        
        dates = [
            p for p in parts
            if 'T' in p and p.endswith('Z') and len(p) == 16
        ]
        
        if len(dates) >= 2:
            pair = (dates[0][:8], dates[1][:8])
            proc_time = dates[2] if len(dates) >= 3 else ""
            
            if pair not in daac_file_map:
                daac_file_map[pair] = []
            daac_file_map[pair].append((proc_time, filename))
            
    for pair in daac_file_map.keys():
        daac_pairs.append(pair)
        
    return sorted(daac_pairs), daac_file_map


def plot_original_2_panel(daac_pairs, new_pairs, frame_id):
    """Generates the 2-panel comparison and saves to PNG."""
    new_dates = set(
        datetime.strptime(d, "%Y%m%d") for p in new_pairs for d in p
    )
    daac_dates = set(
        datetime.strptime(d, "%Y%m%d") for p in daac_pairs for d in p
    )
    
    min_overlap = max(min(new_dates), min(daac_dates))
    max_overlap = min(max(new_dates), max(daac_dates))

    daac_only_ovlp = {
        d for d in (daac_dates - new_dates)
        if min_overlap <= d <= max_overlap
    }
    new_only_ovlp = {
        d for d in (new_dates - daac_dates)
        if min_overlap <= d <= max_overlap
    }

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    title = f"Acquisition Timeline & Verification - Frame: {frame_id}"
    fig.suptitle(title, fontsize=16, fontweight='bold')

    pair_counts = Counter(new_pairs)
    unique_new_pairs = sorted(list(pair_counts.keys()))
    duplicates_found = []

    gap_idx = 1
    for i, (start, end) in enumerate(unique_new_pairs):
        s_date = datetime.strptime(start, "%Y%m%d")
        e_date = datetime.strptime(end, "%Y%m%d")
        
        ax1.hlines(y=i, xmin=s_date, xmax=e_date, color='blue', linewidth=4)
        
        count = pair_counts[(start, end)]
        if count > 1:
            ax1.plot(
                e_date, i, marker='*', color='gold', 
                markeredgecolor='darkorange', markersize=10, zorder=5
            )
            dup_label = (
                f"{start[:4]}-{start[4:6]}-{start[6:]} to "
                f"{end[:4]}-{end[4:6]}-{end[6:]} ({count}x)"
            )
            duplicates_found.append(dup_label)
        
        if i < len(unique_new_pairs) - 1:
            next_start = datetime.strptime(unique_new_pairs[i+1][0], "%Y%m%d")
            if next_start > e_date:
                gap_lbl = (
                    f"Gap {gap_idx}: {e_date.strftime('%Y-%m-%d')} "
                    f"to {next_start.strftime('%Y-%m-%d')}"
                )
                ax1.axvspan(
                    e_date, next_start, color='red', alpha=0.2, label=gap_lbl
                )
                
                mid_date = e_date + (next_start - e_date) / 2
                ax1.annotate(
                    str(gap_idx), (mid_date, i), xytext=(0, 10),
                    textcoords='offset points', color='red',
                    fontsize=10, fontweight='bold', ha='center'
                )
                gap_idx += 1
    
    ax1.set_title("SDS/forward mode products")
    ax1.set_yticks([])
    ax1.grid(axis='x', linestyle='--', alpha=0.7)
    
    if duplicates_found:
        ax1.plot(
            [], [], marker='*', color='gold', markeredgecolor='darkorange', 
            linestyle='None', markersize=10, label='Duplicate Products'
        )
        
    handles, labels = ax1.get_legend_handles_labels()
    if handles:
        ax1.legend(
            handles, labels, loc='center left', 
            bbox_to_anchor=(1.02, 0.5), title="Timeline Flags"
        )

    if duplicates_found:
        dup_str = "\n".join(duplicates_found[:8])
        if len(duplicates_found) > 8:
            dup_str += f"\n... (+{len(duplicates_found)-8} more)"
        
        ax1.annotate(
            f"Duplicate Products:\n{dup_str}", 
            xy=(0.02, 0.95), xycoords='axes fraction',
            fontsize=8, va='top', ha='left',
            bbox=dict(
                boxstyle="round,pad=0.4", fc="lightyellow", 
                ec="gold", alpha=0.8
            )
        )

    all_spans = []
    for p in daac_pairs + new_pairs:
        s_dt = datetime.strptime(p[0], "%Y%m%d")
        e_dt = datetime.strptime(p[1], "%Y%m%d")
        all_spans.append((s_dt, e_dt))
    
    all_spans.sort(key=lambda x: x[0])
    
    combined_gaps = []
    if all_spans:
        curr_max = all_spans[0][1]
        for s_dt, e_dt in all_spans[1:]:
            if s_dt > curr_max:
                combined_gaps.append((curr_max, s_dt))
                curr_max = e_dt
            else:
                curr_max = max(curr_max, e_dt)

    cg_added = False
    for g_start, g_end in combined_gaps:
        lbl = "Network Connectivity Gap" if not cg_added else ""
        ax2.axvspan(
            g_start, g_end, color='red', alpha=0.3, hatch='//', label=lbl
        )
        mid_date = g_start + (g_end - g_start) / 2
        ax2.annotate(
            f"Connectivity Gap:\n{g_start.strftime('%Y-%m-%d')}\n"
            f"to {g_end.strftime('%Y-%m-%d')}",
            (mid_date, 0.5), ha='center', va='center',
            fontsize=8, color='darkred', fontweight='bold',
            bbox=dict(
                boxstyle="round,pad=0.2", fc="mistyrose", 
                ec="red", alpha=0.9
            )
        )
        cg_added = True

    ax2.scatter(
        list(daac_dates), [0]*len(daac_dates),
        color='red', s=100, label="DAAC/Historic products", zorder=3
    )
    ax2.scatter(
        list(new_dates), [1]*len(new_dates),
        color='blue', s=100, label="SDS/forward mode products", zorder=3
    )

    target_str_dates = ['20250429', '20250430', '20250501']
    target_dates = {
        datetime.strptime(d, "%Y%m%d") for d in target_str_dates
    }
    degraded_sds = new_dates.intersection(target_dates)

    if degraded_sds:
        ax2.scatter(
            list(degraded_sds), [1] * len(degraded_sds),
            color='purple', marker='D', s=100, zorder=5,
            label=("Products with potential degraded burst "
                   "synchronization")
        )
    
    if daac_only_ovlp:
        ax2.scatter(
            list(daac_only_ovlp), [0]*len(daac_only_ovlp),
            color='orange', marker='X', s=150, zorder=4,
            label="DAAC/Historic Only (Missing in SDS/forward)"
        )
        for i, d in enumerate(sorted(daac_only_ovlp)):
            ax2.annotate(
                d.strftime("%Y-%m-%d"), xy=(d, 0),
                xytext=(-15, -25 - (i % 3) * 15), textcoords='offset points',
                color='darkorange', fontweight='bold', ha='right', va='top',
                arrowprops=dict(
                    arrowstyle='-', color='darkorange', alpha=0.7
                )
            )
            
    if new_only_ovlp:
        ax2.scatter(
            list(new_only_ovlp), [1]*len(new_only_ovlp),
            color='green', marker='X', s=150, zorder=4,
            label="SDS/forward Only (Missing in DAAC)"
        )
        for i, d in enumerate(sorted(new_only_ovlp)):
            ax2.annotate(
                d.strftime("%Y-%m-%d"), xy=(d, 1),
                xytext=(-15, 25 + (i % 3) * 15), textcoords='offset points',
                color='darkgreen', fontweight='bold', ha='right', va='bottom',
                arrowprops=dict(
                    arrowstyle='-', color='darkgreen', alpha=0.7
                )
            )

    all_dates = new_dates.union(daac_dates)
    for d in all_dates:
        ax2.vlines(
            x=d, ymin=-0.5, ymax=1.5,
            color='gray', linestyle=':', alpha=0.5, zorder=1
        )

    display_min = min(new_dates)
    display_max = max(new_dates)
    
    for g_start, g_end in combined_gaps:
        if g_end >= display_min and g_start <= display_max:
            display_min = min(display_min, g_start)
            display_max = max(display_max, g_end)

    pad = timedelta(days=max(5, (display_max - display_min).days // 20))
    ax2.set_xlim(display_min - pad, display_max + pad)

    ax2.set_title("Mismatch Check (Colored 'X' denotes missing overlap)")
    ax2.set_yticks([0, 1])
    ax2.set_yticklabels(["DAAC/Historic List", "SDS/forward sample products"])
    ax2.set_ylim(-0.5, 1.5)
    
    ax2.legend(loc='center left', bbox_to_anchor=(1.02, 0.5))
    ax2.grid(axis='x', linestyle='--', alpha=0.7)

    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    fig.autofmt_xdate(rotation=45, ha='right')

    plt.tight_layout()
    plt.subplots_adjust(right=0.72)

    out_file = f"acquisition_timeline_2panel_{frame_id}.png"
    plt.savefig(out_file, dpi=300, bbox_inches='tight')
    print(f"2-Panel Plot saved to: {out_file}")
    plt.close(fig)


def plot_reference_vs_sec_scatter(daac_pairs, s3_pairs, frame_id):
    """Plots Ref vs Sec dates and explicitly labels unshared dates."""
    fig, ax = plt.subplots(figsize=(10, 10))
    title = f"Ministacking Structure - Frame: {frame_id}"
    ax.set_title(title, fontsize=16, fontweight='bold')

    daac_ref = [datetime.strptime(p[0], "%Y%m%d") for p in daac_pairs]
    daac_sec = [datetime.strptime(p[1], "%Y%m%d") for p in daac_pairs]

    s3_ref = [datetime.strptime(p[0], "%Y%m%d") for p in s3_pairs]
    s3_sec = [datetime.strptime(p[1], "%Y%m%d") for p in s3_pairs]

    ax.scatter(
        daac_ref, daac_sec,
        marker='D', s=100, facecolors='none', edgecolors='red',
        linewidths=1.2, label=f"DAAC/Historic ({len(daac_pairs)})", zorder=1
    )

    s3_label = f"SDS/forward mode products ({len(s3_pairs)})"
    ax.scatter(
        s3_ref, s3_sec,
        marker='o', s=30, color='blue',
        label=s3_label, zorder=2
    )

    daac_unique = set(d for p in daac_pairs for d in p)
    s3_unique = set(d for p in s3_pairs for d in p)
    
    min_overlap = max(min(daac_unique), min(s3_unique))
    max_overlap = min(max(daac_unique), max(s3_unique))
    
    daac_only = {
        d for d in (daac_unique - s3_unique) 
        if min_overlap <= d <= max_overlap
    }
    s3_only = {
        d for d in (s3_unique - daac_unique) 
        if min_overlap <= d <= max_overlap
    }
    
    if daac_only or s3_only:
        mm_ref, mm_sec = [], []
        annotations = {}
        marked_dates = set()
        
        for r, s in daac_pairs + s3_pairs:
            mark_this_pair = False
            labels = []
            
            if r in daac_only and r not in marked_dates:
                labels.append(f"{r[:4]}-{r[4:6]}-{r[6:]} (DAAC/Historic only)")
                marked_dates.add(r)
                mark_this_pair = True
            elif r in s3_only and r not in marked_dates:
                labels.append(f"{r[:4]}-{r[4:6]}-{r[6:]} (SDS/forward only)")
                marked_dates.add(r)
                mark_this_pair = True
                
            if s in daac_only and s not in marked_dates:
                labels.append(f"{s[:4]}-{s[4:6]}-{s[6:]} (DAAC/Historic only)")
                marked_dates.add(s)
                mark_this_pair = True
            elif s in s3_only and s not in marked_dates:
                labels.append(f"{s[:4]}-{s[4:6]}-{s[6:]} (SDS/forward only)")
                marked_dates.add(s)
                mark_this_pair = True
            
            if mark_this_pair:
                ref_dt = datetime.strptime(r, "%Y%m%d")
                sec_dt = datetime.strptime(s, "%Y%m%d")
                mm_ref.append(ref_dt)
                mm_sec.append(sec_dt)
                
                if (ref_dt, sec_dt) in annotations:
                    annotations[(ref_dt, sec_dt)] += "\n" + "\n".join(labels)
                else:
                    annotations[(ref_dt, sec_dt)] = "\n".join(labels)
                
        if mm_ref:
            ax.scatter(
                mm_ref, mm_sec,
                marker='x', s=80, color='orange', linewidths=2,
                label=f"Unshared date (in overlap window)", zorder=3
            )
            
            stagger_positions = [
                (-30, 20, 'right', 'bottom'),  
                (20, -30, 'left', 'top'),      
                (-50, 40, 'right', 'bottom'),  
                (40, -50, 'left', 'top')       
            ]
            
            off_idx = 0
            for (x, y), text in annotations.items():
                dx, dy, align_h, align_v = stagger_positions[
                    off_idx % len(stagger_positions)
                ]
                ax.annotate(
                    text, (x, y), xytext=(dx, dy), textcoords='offset points',
                    fontsize=8, color='darkorange', fontweight='bold',
                    ha=align_h, va=align_v,
                    arrowprops=dict(
                        arrowstyle='-', color='darkorange', alpha=0.7
                    ),
                    bbox=dict(
                        boxstyle="round,pad=0.2", fc="white", 
                        alpha=0.7, ec="none"
                    )
                )
                off_idx += 1

    ax.margins(0.15)
    
    ax.set_xlabel("Reference (Start Date)", fontsize=12)
    ax.set_ylabel("Secondary (End Date)", fontsize=12)
    
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    ax.yaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    fig.autofmt_xdate(rotation=45, ha='right')

    ax.grid(True, linestyle='-', alpha=0.7)
    ax.legend(loc='upper left', fontsize=12)

    plt.tight_layout()
    out_file = f"disp_ref_vs_sec_{frame_id}.png"
    plt.savefig(out_file, dpi=300, bbox_inches='tight')
    print(f"Ref vs Sec Scatter Plot saved to: {out_file}")
    plt.close(fig)


if __name__ == "__main__":
    args = parse_args()
    
    # 1. First attempt to load SDS data locally
    sds_data, sds_file_map, sds_target_dir = fetch_sds_local_pairs(args.sds_local_dir, args.frame)
    
    # 2. Check if local data was found; if not, fallback to S3
    if sds_data:
        print(f"Found {len(sds_data)} SDS products locally in '{sds_target_dir}'. Prioritizing local data.")
    else:
        print(f"No local SDS products found in '{args.sds_local_dir}'. Falling back to S3 search...")
        sds_data, sds_file_map = fetch_s3_pairs(args.frame, args.s3_path)
        
    daac_data, daac_file_map = parse_daac_directory(args.daac_dir, args.frame)
    
    if not sds_data and not daac_data:
        print("No pairs found from either source. Exiting.")
        sys.exit(1)

    # Generate both independent plots
    if sds_data and daac_data:
        plot_original_2_panel(daac_data, sds_data, args.frame)
        plot_reference_vs_sec_scatter(daac_data, sds_data, args.frame)
        
    daac_unique_dates = set(d for pair in daac_data for d in pair)
    sds_unique_dates = set(d for pair in sds_data for d in pair)
    
    daac_min = min(daac_unique_dates) if daac_unique_dates else "N/A"
    daac_max = max(daac_unique_dates) if daac_unique_dates else "N/A"
    sds_min = min(sds_unique_dates) if sds_unique_dates else "N/A"
    sds_max = max(sds_unique_dates) if sds_unique_dates else "N/A"
    
    print(
        f"\nNumber of DAAC/Historic mode products: {len(daac_data)} "
        f"({len(daac_unique_dates)} unique acquisition dates. "
        f"Earliest: {daac_min}, Most recent: {daac_max})"
    )
    print(
        f"Number of unique SDS/forward mode products: {len(sds_file_map)} "
        f"({len(sds_unique_dates)} unique acquisition dates. "
        f"Earliest: {sds_min}, Most recent: {sds_max})"
    )
    
    target_str_dates = {'20250429', '20250430', '20250501'}
    flagged = []
    for pair, filenames in sds_file_map.items():
        if pair[0] in target_str_dates or pair[1] in target_str_dates:
            for fname in filenames:
                flagged.append((pair, fname))
                
    if flagged:
        print("\nProducts with potential degraded burst synchronization:")
        for pair, fname in flagged:
            print(f"  Acquisitions: {pair[0]}_{pair[1]} | File: {fname}")
            
    duplicates_to_print = []
    for pair, filenames in sds_file_map.items():
        if len(filenames) > 1:
            duplicates_to_print.append((pair, filenames))
            
    if duplicates_to_print:
        print("\nDuplicate SDS/forward mode products found:")
        for pair, filenames in duplicates_to_print:
            print(f"  Acquisitions: {pair[0]}_{pair[1]}")
            for fname in filenames:
                print(f"    - {fname}")

    daac_duplicates_to_print = []
    for pair, files in daac_file_map.items():
        if len(files) > 1:
            sorted_files = sorted(files, key=lambda x: x[0])
            daac_duplicates_to_print.append((pair, sorted_files))
            
    if daac_duplicates_to_print:
        print("\nDuplicate DAAC/Historic mode products found:")
        for pair, sorted_files in daac_duplicates_to_print:
            print(f"  Acquisitions: {pair[0]}_{pair[1]}")
            for proc, fname in sorted_files:
                note = (
                    " (Selected/Most Recent)" 
                    if fname == sorted_files[-1][1] else ""
                )
                print(f"    - {fname}{note}")

    # --- MINIMUM SPANNING NETWORK SYMLINKS ---
    print(f"\nCreating minimum spanning network in '{args.sds_daac_combo}'...")
    os.makedirs(args.sds_daac_combo, exist_ok=True)
    
    pivot_date = min(sds_unique_dates) if sds_unique_dates else None
    
    if pivot_date:
        selected_daac_pairs = [p for p in daac_data if p[1] <= pivot_date]
        selected_sds_pairs = [p for p in sds_data if p[0] >= pivot_date]
    else:
        selected_daac_pairs = daac_data
        selected_sds_pairs = sds_data

    linked_daac = 0
    for pair in selected_daac_pairs:
        sorted_files = sorted(daac_file_map[pair], key=lambda x: x[0])
        most_recent_filename = sorted_files[-1][1]
        
        src_path = os.path.abspath(
            os.path.join(args.daac_dir, most_recent_filename)
        )
        dest_path = os.path.join(args.sds_daac_combo, most_recent_filename)
        
        if not os.path.lexists(dest_path):
            try:
                os.symlink(src_path, dest_path)
            except Exception as e:
                print(f"Error linking DAAC file {most_recent_filename}: {e}")
                continue
        linked_daac += 1

    linked_sds = 0
    for pair in selected_sds_pairs:
        sorted_sds_files = sorted(sds_file_map[pair])
        most_recent_filename = sorted_sds_files[-1]
        
        # Now safely points to either the local SDS cache or the designated sds_target_dir
        src_path = os.path.abspath(
            os.path.join(sds_target_dir, most_recent_filename)
        )
        dest_path = os.path.join(args.sds_daac_combo, most_recent_filename)
        
        if not os.path.lexists(dest_path):
            try:
                os.symlink(src_path, dest_path)
            except Exception as e:
                print(f"Error linking SDS file {most_recent_filename}: {e}")
                continue
        linked_sds += 1

    print(f"  - DAAC/Historic products linked: {linked_daac}")
    print(f"  - SDS/forward products linked: {linked_sds}")
    
    sys.exit(0)
