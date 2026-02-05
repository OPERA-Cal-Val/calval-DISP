#!/usr/bin/env python
"""
Access and compare metadata for OPERA L3 DISP-S1 NetCDF files.
"""

import argparse
import json
import sys

import xarray as xr
import yaml
from deepdiff import DeepDiff


def get_metadata(file_path):
    """
    Open a DISP-S1 file and return parsed metadata components.
    """
    try:
        ds_meta = xr.open_dataset(file_path, group='metadata')

        # get and parse params
        algo_params_raw = ds_meta.algorithm_parameters_yaml.item()
        runconfig_raw = ds_meta.pge_runconfig.item()
        algo_dict = yaml.safe_load(algo_params_raw)
        run_dict = yaml.safe_load(runconfig_raw)

        # get versions
        sas_version = ds_meta.disp_s1_software_version.item()
        dolphin_version = ds_meta.dolphin_software_version.item()
        compass_version = ds_meta.source_data_software_COMPASS_version.item()
        isce3_version = ds_meta.source_data_software_ISCE3_version.item()
        s1read_version = ds_meta.source_data_software_s1_reader_version.item()

        return {
            "algo": algo_dict,
            "runconfig": run_dict,
            "sas_version": sas_version,
            "dolphin_version": dolphin_version,
            "compass_version": compass_version,
            "isce3_version": isce3_version,
            "s1read_version": s1read_version
        }
    except Exception as err:
        print(f"Error processing {file_path}: {err}")
        sys.exit(1)


def clean_opera_path(path_str):
    """
    Truncates an OPERA file path to remove the processing datetime.
    """
    if not isinstance(path_str, str) or ".h5" not in path_str:
        return path_str
    
    parts = path_str.split('_')
    # OPERA pattern: ..._PROCTIME_POL_VERSION.h5
    # Strip the last 3 elements (time, pol, version)
    if len(parts) > 3:
        return "_".join(parts[:-3])
    return path_str


def cslc_compare_func(obj, path):
    """
    DeepDiff callback for older versions/specific configs.
    Takes 2 arguments: the object being evaluated and its path string.
    Returns True if the object should be EXCLUDED from the diff.
    """
    # If the path contains our target list
    if "cslc_file_list" in path:
        # We return True to EXCLUDE it if the 'cleaned' version 
        # matches the reference. However, DeepDiff's exclude_obj_callback 
        # usually checks if the current obj should be ignored entirely.
        # To ignore "changes" but keep "additions/deletions", 
        # we check the path and return True.
        return True
    return False


def main():
    """Main execution entry point."""
    parser = argparse.ArgumentParser(
        description="Access and compare OPERA L3 DISP-S1 metadata."
    )
    parser.add_argument(
        "file_path",
        help="Path to the primary OPERA DISP-S1 NetCDF file."
    )
    parser.add_argument(
        "--compare",
        help="Path to a second file to compare against the first.",
        default=None
    )

    args = parser.parse_args()
    data1 = get_metadata(args.file_path)

    if args.compare:
        data2 = get_metadata(args.compare)

        print(f"Comparing File 1: {args.file_path}")
        print(f"Against File 2: {args.compare}\n")
        print("=" * 60)

        # 1. Compare Software Version
        checks = [
            ('sas_version', 'SAS version'),
            ('dolphin_version', 'Dolphin software version'),
            ('compass_version', 'COMPASS software version'),
            ('isce3_version', 'ISCE3 software version'),
            ('s1read_version', 'S1 reader software version'),
        ]

        for key, label in checks:
            val1 = data1.get(key)
            val2 = data2.get(key)

            if val1 == val2:
                print(f"✅ {label}: Match ({val1})")
            else:
                print(f"❌ {label}: Mismatch {val1} vs {val2}")

        # 2. Compare Dictionaries (Algo and Runconfig)
        for key in ['algo', 'runconfig']:
            print(f"\n--- {key.upper()} DIFFERENCES ---")
            
            # For your specific case of ignoring the processing time 
            # in a list of 1000+ files, we'll exclude the cslc_file_list 
            # from the standard diff and check it separately if you wish, 
            # or simply exclude the objects via callback.
            diff = DeepDiff(
                data1[key], 
                data2[key], 
                ignore_order=True,
                exclude_obj_callback=cslc_compare_func
            )

            if not diff:
                print(f"No significant differences found in {key}.")
            else:
                print(diff.pretty())
    else:
        print("### ALGORITHM PARAMETERS ###")
        print(json.dumps(data1['algo'], indent=4))
        print("\n### RUNCONFIG ###")
        print(json.dumps(data1['runconfig'], indent=4))
        print(f"\ndolphin software version: {data1['dolphin_version']}")
        print(f"\nsas version: {data1['sas_version']}")


if __name__ == "__main__":
    main()

