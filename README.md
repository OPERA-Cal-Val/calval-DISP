# OPERA Surface Displacement (DISP) validation tools
Repository of tools intended to address the Cal/Val validation plan for the Level 3 OPERA DISP product suite. Namely, verify that products meet the established OPERA L3 DISP uncertainty requirements.

## Contents

1. [Installation](#installation)
2. [General Procedures](#general-procedures)
    - [InSAR LOS velocity estimation](#insar-los-velocity-estimation)
    - [Validation data LOS velocity estimation](#validation-data-los-velocity-estimation)
    - [Statistical analyses](#statistical-analyses)
3. [Validation Approaches](#validation-approaches)
    - [Validation approach 1 (VA1)](#validation-approach-1-va1)
    - [Validation approach 2 (VA2)](#validation-approach-2-va2)
    - [Validation approach 3 (VA3) TBD](#validation-approach-3-va3-tbd)
4. [Running Validation Workflow Guide](#running-validation-workflow-guide)
    - [Step 0: Download GNSS measurements](#step-0-download-gnss-measurements)
    - [Step 1: Download DISP-S1 products, static layers, and geometry files](#step-1-download-disp-s1-products-static-layers-and-geometry-files)
    - [Step 2: Prepare timeseries and velocity estimate](#step-2-prepare-timeseries-and-velocity-estimate)
    - [Step 3: Generating non-linear displacement score map](#step-3-generating-non-linear-displacement-score-map)
    - [Step 4: Perform VA1 and VA2 Cal/Val analysis](#step-4-perform-va1-and-va2-calval-analysis)
    - [Step 5: OPTIONAL Generate comprehensive GNSS vs InSAR comparison plots](#step-5-optional-generate-comprehensive-gnss-vs-insar-comparison-plots)
5. [Validation Report Tool](#validation-report-tool)
6. [OPTIONAL Compare DISP Layers Between Processing Versions](#optional-compare-disp-layers-between-processing-versions)
7. [OPTIONAL dolphin Software Workflow](#optional-dolphin-software-workflow)
8. [Legacy Tools](#legacy-tools)
9. [References](#references)
10. [Contributors](#contributors)

## Installation

A detailed installation on how to build a working environment can be found __[here](https://github.com/OPERA-Cal-Val/calval-DISP/blob/main/docs/installation.md/)__.
<br />
<br />
<br />
<br />

## General Procedures

Summarized below are the general procedures to implement our validation approaches (outlined in the next section) including the two primary statistical analyses that are performed:

### InSAR LOS velocity estimation

We estimate the long term line-of-sight rates and corresponding uncertainties by fitting a time-series function which includes a linear rate parameter across the displacement history, derived from the reconstructed displacement time-series using consecutive DISP products in time. Where needed, we may include additional complexities (e.g. Heaviside step, exponential or logarithmic function for coseismic and postseismic events/processes).

### Validation data LOS velocity estimation

We project the independent Validation displacement data (e.g GNSS time-series) in the radar look direction and then estimate the velocity over the same period of time and with the same time-series models as from step 1. The GNSS time-series data will have its processing and noise corrections applied. As an example, the GNSS time-series data obtained from the University of Nevada, Reno (UNR) is corrected for instrument changes, ionospheric and tropospheric propagation delays, as well as solid earth tides [Blewitt et al., 2018].

### Statistical analyses

We mainly implement similar analyses as exemplified to assess if the OPERA L3 DISP product requirements are achieved.
  - a)	Variogram analysis for the line-of-sight velocity residuals will be performed. Using the structure function, values at length scales 0.1 km < L < 50 km will be extracted and then compared to the OPERA L3 DISP requirements for validation.
  - b)	We calculate the mean and the standard deviation of the velocity residuals (e.g. perform scatter plot/RMSE analysis) and then assess if the OPERA L3 DISP requirements are satisfied.

For all the validation analyses, we report the accuracy that can be achieved with and without applying the optional correction layers embedded within the L3 DISP products, including the ionospheric propagation delays (e.g., [Fattahi et al., 2017]), tropospheric propagation delays (e.g., [Jolivet et al., 2011] and [Doin et al., 2009]), and solid earth tide effects (e.g. [Agnew 2012]).
<br />
<br />
<br />
<br />

## Validation Approaches

Three validation approaches will be performed separately depending on the Validation dataset are used to assess the DISP product requirements. Each validation approach is dubbed in this document as Validation Approach n (VAn) for n = 1, 2, 3, i.e. VA1, VA2, and VA3. The validation approaches are primarily differentiated based on the surface displacement occurring at the validation site (i.e. deforming or negligibly deforming) and/or the instrument used to validate the InSAR measurement (i.e. GNSS, corner reflectors, InSAR). Please refer to the __[OPERA validation plan](https://d2pn8kiwq2w21t.cloudfront.net/documents/OPERA_Validation_Plan.pdf) for a more exhaustive overview.

### Validation approach 1 (VA1)

Perform direct comparison of GNSS and InSAR is performed over deforming regions with dense GNSS networks including Central Valley in California, Hawaii, New York City, Hampton Roads in Virginia, Seattle, and Houston/Galveston. These validation sites are shown in Figure 1 and detailed in Table 1. Temporally decorrelated pixels on the InSAR data (i.e. temporal coherence ≤ 0.5) are removed and a water mask is applied during validation. An initial inspection of the GNSS time-series and the study area is done to select a stable reference area for all the InSAR data (i.e. reference point noise analysis). We perform root-mean-square error (RMSE) analysis (i.e. statistical analysis 3b above) to assess if we meet our DISP product requirements.

<p align="center">
  <img width="90%" height="90%" src="https://github.com/OPERA-Cal-Val/calval-DISP/assets/13227405/a32159f3-ae35-4a02-b174-83247a23cd60">
</p>
<B>Figure 1.</B> Proposed Validation displacement sites from Table 4.3 displayed geographically. Red outlines are regions which have dense GNSS networks and are considered deforming. Orange outline (Oklahoma) is considered a negligibly deforming region with dense GNSS stations. Gray square in California is the target site of the CRs to be deployed by OPERA.
<br />
<br />


<p align="center">
  <img width="90%" height="90%" src="https://github.com/OPERA-Cal-Val/calval-DISP/assets/13227405/5bf73799-ff2f-451a-8408-bdb6c3636d31">
</p>
<B>Table 1.</B> Proposed Validation displacement sites chosen to represent a diversity of deformation processes in Northern America over various climate and terrain locations.
<br />
<br />

### Validation approach 2 (VA2)

InSAR residual analysis is conducted over a variety of negligibly deforming regions or in areas that are located outside the deformation field within the SAR footprint, including the Mojave Desert in California and Oklahoma. For negligibly deforming regions, we target those that have maintained coherence (i.e. temporal coherence > 0.5) PS and DS pixels within a frame. For a scene that includes a deformation field, we will mask first the deformation area as well as those that have temporal coherence ≤ 0.5 prior to residual analysis.

The InSAR displacement rate residuals will be computed from a dense grid of randomly selected pairs of pixels to evaluate the accuracy as a function of distance. Since the selection of pixels is unbiased (or random), the inherent accuracy is evaluated not only on PS pixels but on a combined performance of both PS and DS within the SAR footprint. Note that this VA2 strategy supports and complements VA1 (e.g. length scales, L < 0.1 km ) and also removes the limitation in terms of the number of GNSS stations that can be used to validate the product requirements. We perform a variogram analysis as described above (i.e. statistical analysis 3a) to assess if we meet the DISP requirements (i.e. 3 mm per yr or better for Sentinel-1 and 5 mm per yr or better for NISAR, over length scales of 0.1 km < L < 50 km as reported in Table 2).

The North America Sentinel-1 12 day coherence study previously done for NISAR (Figure 4.10) will be used as a reference guide to ensure that the selected negligibly deforming validation sites meet the temporal coherence requirement.

<p align="center">
  <img width="90%" height="90%" src="https://github.com/OPERA-Cal-Val/calval-DISP/assets/13227405/27ec413a-73c0-4fe7-aa67-815b1bfc46ca">
</p>
<B>Table 2.</B> Summary of OPERA DISP accuracy requirements.
<br />
<br />

### Validation approach 3 (VA3) TBD

In addition to the two validation approaches described above, we are also performing a validation experiment using the corner reflectors that we will deploy across the creeping segment of the San Andreas Fault. The precise height and locations of the CRs are obtained from their co-located GNSS stations or by conducting geodetic surveys annually or semi-annually. We will take the 3-component CR position values and project them to the satellite’s line-of-sight direction to estimate the CR LOS velocities. We will then calculate the residuals by directly comparing the CR LOS velocities and InSAR LOS displacement rates. To test if the requirements are achieved for NISAR and Sentinel-1, we will implement statistical analysis 3b as described above. Additional rigorous statistical modeling and tests will be considered for additional verification of the requirements. NOTE: this approach has not yet been implemented in this repository.
<br />
<br />
<br />
<br />

## Running Validation Workflow Guide

Outlined below are the components of the validation workflow and how to run and manage each step.

We outline a sample run for the Central Valley, California case study for descending Sentinel-1 track 042, which maps roughly to an OPERA frame ID of 11116.

For a list of all Cal/Val sites and corresponding basic input parameters, refer to [DISP-S1_CalVal_sites.csv](https://github.com/OPERA-Cal-Val/calval-DISP/blob/main/validation_data/DISP-S1_CalVal_sites.csv).

### Step 0: Download GNSS measurements

Using the [run0_gnss_download_screen.py](https://github.com/OPERA-Cal-Val/calval-DISP/blob/main/run0_gnss_download_screen.py) script, download prerequisite GNSS measurements needed for Cal/Val.
```bash
# Args:
# --frameID    OPERA frame number
# --startDate  Start date (default: 20160701)
# --endDate    End date (default: 20240930)
# --gnss_completeness_threshold    Ratio of data timespan with valid GNSS epochs (default: 0.8)
# --gnss_thr_eq                    Threshold of earthquake magnitude for removing GNSS stations (default: 11)
# --gnss-source                    GNSS source repo (default: UNR)

cd /path/to/work/folder
run0_gnss_download_screen.py \
      --frameID 11116
```
<br />
<br />

### Step 1: Download DISP-S1 products, static layers, and geometry files

Using the [run1_download_DISP_S1_Static.py](https://github.com/OPERA-Cal-Val/calval-DISP/blob/main/run1_download_DISP_S1_Static.py) script, download DISP-S1 products and prerequisite static layers and geometry files needed for Cal/Val.

By default, the script downloads all available dates. You can select a specific date range using the `--startDate` and `--endDate` arguments.
```bash
# Args:
# --frameID    OPERA frame number
# --version    OPERA dataset version, defaults to production version 1.0 (default: 1.0)
# --staticDir  Folder for static layers/metadata
# --geomDir    Folder for geometry files
# --dispDir    Folder for data
# --startDate  Start date (optional)
# --endDate    End date (optional)

run1_download_DISP_S1_Static.py \
      --frameID 11116 \
      --staticDir /path/to/work/folder/static_lyrs \
      --geomDir /path/to/work/folder/geometry \
      --dispDir /path/to/work/folder/data
```
<br />
<br />

### Step 2: Prepare timeseries and velocity estimate

Using the [run2_prep_mintpy_opera.py](https://github.com/OPERA-Cal-Val/calval-DISP/blob/main/run2_prep_mintpy_opera.py) script, extract displacement layers and prepare expected timeseries (`timeseries.h5`) and velocity estimate (`velocity.h5`) inputs for Cal/Val analysis. Outputs are extracted in MintPy-software compatible HDF5 file structures.

To access the pre-defined reference points for each of the designated Cal/Val layers, refer again to [DISP-S1_CalVal_sites.csv](https://github.com/OPERA-Cal-Val/calval-DISP/blob/main/validation_data/DISP-S1_CalVal_sites.csv).
```bash
# Args:
# -m   Folder for static layers/metadata
# -u   Folder with data (*.nc for all files)
# -g   Folder for geometry files
# -o   Folder for timeseries output
# --ref-lalo         Spatial reference for timeseries pass as 'Lat Lon'
# --water-mask-file  Water mask file (default: esa_world_cover_2021)
# --dem-file         DEM file (default: glo_30)

run2_prep_mintpy_opera.py \
        -m "/path/to/work/folder/static_lyrs" \
        -u "/path/to/work/folder/data/*.nc" \
        -g "/path/to/work/folder/geometry" \
        -o /path/to/work/folder/mintpy_output \
        --ref-lalo '36.612 -121.064'
```

You can visualize the outputs using MintPy
```bash
## Need help with the arguments: tsview.py -h
tsview.py /path/to/work/folder/mintpy_output/velocity.h5 \
        -m /path/to/work/folder/mintpy_output/recommended_mask90threshold.h5 \
```

Note: 
`recommended_mask90threshold.h5` is based on the time-series of `recommended_mask` layers (i.e. `recommended_mask.h5`). We by default (i.e. `--reliability-threshold 0.9`) pick the top 90% representing the "most reliable pixels in time" after normalizing the `recommended_mask` against the total number of epoch/dataset.
<br />
<br />

### Step 3: Generating non-linear displacement score map

Using the [run3_cal_non_linDisp_score.py](https://github.com/OPERA-Cal-Val/calval-DISP/blob/main/run3_cal_non_linDisp_score.py) script, generate a non-linear displacement score map `nonDispScore.h5`. This is intended to mask the InSAR velocity estimates for Cal/Val to meet the expected VA2 condition of analysis on negligibly deforming regions.
```bash
run3_cal_non_linDisp_score.py
```
<br />
<br />

### Step 4: Perform VA1 and VA2 Cal/Val analysis

The __[run4_DISP-S1_Secular_Requirement_Validation.ipynb](https://github.com/OPERA-Cal-Val/calval-DISP/blob/main/run4_DISP-S1_Secular_Requirement_Validation.ipynb)__ notebook contains a suite of functionality to evaluate the uncertainties associated with the validation requirements as outlined above. The notebook itself was adopted the NISAR team Algorithm Theoretical Basis Document (ATBD) __[notebook for secular motion](https://github.com/nisar-solid/ATBD/blob/main/methods/secular/Secular_Requirement_Validation.ipynb/)__.

End-to-end runs for select sites can be found under the [calval_runs](https://github.com/OPERA-Cal-Val/calval-DISP/tree/main/calval_runs) folder.
```bash
# Args:
# site         DISP-S1 Frame ID prepended with 'F' (e.g. 'F11116')
# gnss_source  GNSS source repo (default: UNR)

# coherence_threshold       Coherence threshold (default: 0.6)
# apply_coh_mask            Apply coherence mask based off threshold above (default: True)
# reliability_threshold     Threshold for reliability based on the counts of valid pixels (default: 0.9)
# apply_recommended_mask    Apply reliability mask based off threshold above (default: True)
# outlier_removal_method    Outlier removal method options 'zscore' or 'modified_zscore' (default: modified_zscore)
# outlier_zscore_threshold  Threshold if 'zscore' is specified above (default: 2.0)
# apply_outlier_removal     Apply outlier removal mask (default: True)

# gnss_completeness_threshold    Ratio of data timespan with valid GNSS epochs (default: 0.8)
# gnss_residual_stdev_threshold  Max threshold standard deviation of residuals to linear GNSS fit (default: 10.)

# thr_var_score          Variability score threshold (default: 0.4)
# apply_nonlinear_mask   Apply nonlinear mask generated from step 3 (default: True)

# step_events_date     Step events for e.g. earthquakes, volcanos specified as YYYYMMDD (default: None)

# copy necessary files
cp /path/to/source/repo/calval-DISP/run4_DISP-S1_Secular_Requirement_Validation.ipynb .
cp -r /u/trappist-r0/ssangha/conda_installation/miniforge/miniforge/calval-DISP/validation_data .

# run notebook using papermill
papermill run4_DISP-S1_Secular_Requirement_Validation.ipynb run4_completed.ipynb \
        -p site "F11116"
```
<br />
<br />

### Step 5: [OPTIONAL] Generate comprehensive GNSS vs InSAR comparison plots

The __[run5_plot_InSARvsGNSS_TS.ipynb](https://github.com/OPERA-Cal-Val/calval-DISP/blob/main/run5_plot_InSARvsGNSS_TS.ipynb)__ notebook expands upon the simple GNSS vs InSAR comparison plots of the notebook from step 4 with more comprehensive details, including histograms of relative differences. NOTE: This step is entirely optional, it is not necessary for validation.
```bash
# copy necessary files
cp /path/to/source/repo/calval-DISP/run5_plot_InSARvsGNSS_TS.ipynb .

# run notebook using papermill
papermill run5_plot_InSARvsGNSS_TS.ipynb run5_completed.ipynb \
        -p site "F11116"
```
<br />
<br />
<br />
<br />

## Validation Report Tool

The automated generation of the formal validation report used to summarize our Cal/Val activities is supported by our [validation report script](https://github.com/OPERA-Cal-Val/calval-DISP/blob/main/validation_report/prep_disp_validation_report.py).

To use it, run the Validation Workflow as outlined above through to completion for at least one track `FXXXX`.

Your working directory should contain the Cal/Val results of all tracks you wish to summarize. E.g. `/path/to/working/directory` contains `F08622/results`, `F11115/results`, etc.

To run this tool, simply launch it as so:
```bash
/path/to/source/repo/calval-DISP/validation_report/prep_disp_validation_report.py -p /path/to/working/directory --fname calval_report.pdf
```

Inspect the output `calval_report.pdf` to see the validation report summary.
<br />
<br />
<br />
<br />

## [OPTIONAL] Compare DISP Layers Between Processing Versions

This optional QA step uses the `DISP_comparison_suite.ipynb` notebook to
visualize multiple MintPy layers (timeseries, ionosphere, short wavelength
corrections, velocity, etc.) for two processing versions of the **same frame**.
It is intended to be run **after completing Step 2 (“Prepare timeseries and
velocity estimate”)** for at least two versions (for example, `OPS_data` and
`SDS_data`). The notebook relies on the MintPy outputs generated by Step 2 (one
directory per `Fxxxxx` frame under each processing version).

### Configure papermill parameters

Populate the following parameters when launching the suite with papermill:

| Parameter | Description |
| --- | --- |
| `frame_num` | Frame number (e.g., `"08882"`). |
| `parent_dir` | Root directory containing version subfolders produced by Step 2. |
| `version_num` | Name of the “primary” processing version (e.g., `"OPS_data"`). |
| `version_num_to_comp` | Name of the comparison version (e.g., `"SDS_data"`). |
| `threshold_density` | Density threshold used for Section 3 (set to `0.9` by default). |
| `list_lat_lon` | Either a list of `(lat, lon)` tuples, a WKT path, or an integer specifying how many random points to seed. |
| `list_lat_lon_wkt_output` | Optional output path for the canonical list of selected points; defaults to `<parent_dir>/figures/Fxxxxx_selected_points.wkt`. |
| `ref_lalo` | Specify spatial reference point as 'Lat Lon' and re-reference all inputs to this point (defaults to native reference point of inputs). |
| `pixel_radius`, `random_seed`, `timeseries_filename`, `bins`, etc. | Advanced options used by the random-timeseries workflow (pre-populated in the parameter cell). |

> **Tip:** If you pass an integer to `list_lat_lon`, the notebook will pick that
> many valid pixels shared by both stacks and persist them back to WKT so the
> same set can be reused later.

### Execute via papermill

Run papermill from the repository root (or wherever the notebook lives), passing
the parameter overrides as needed:

```bash
papermill \
  DISP_comparison_suite.ipynb \
  run_DISP_comparison_suite.ipynb \
  -p frame_num "08882" \
  -p parent_dir "/path/to/calval_output" \
  -p version_num "OPS_data" \
  -p version_num_to_comp "SDS_data" \
  -p list_lat_lon "/path/to/seed_points.wkt" \
  -p threshold_density 0.9
```

The executed notebook (`run_DISP_comparison_suite.ipynb`) serves as the QA log
for the frame/version comparison.

### Outputs and directory layout

The notebook writes all figures beneath `<parent_dir>/<figures_dir>/`. Each
layer-specific section places its PNGs in a dedicated subdirectory, for example:

```
<figures_dir>/
  displacement/
    compare_detailed_OPS_data_SDS_data_01.png
    displacement_fallback_height.png
    displacement_fallback_points.wkt
  short_wavelength_displacement/
  ionospheric_delay/
  ...
```

- “Compare detailed” PNGs contain the paired timeseries plots for each selected
  pixel.
- Every layer records a fallback WKT if a requested coordinate had to be
  replaced inside that section (e.g., `displacement_fallback_points.wkt`).
- Whenever fallback points are used, a companion height map
  (`*_fallback_height.png`) is generated with those pixels highlighted in red so
  you can visually confirm where new points were sampled.
- The canonical point list selected in Section 4.0 is saved as
  `Fxxxxx_selected_points.wkt` in the `figures` directory (or the path supplied
  via `list_lat_lon_wkt_output`).

### What the notebook does

1. **Section 3:** Recreates the timeseries-density histograms (both thresholded
   and unthresholded) for the two versions, mirroring the Step 4 QA workflow.
2. **Section 4:** Displays epoch-based MintPy layers side-by-side (recommended
   mask, displacement, short-wavelength corrections, ionosphere, solid Earth
   tide, connected components, temporal coherence, estimated phase quality,
   water mask, inversion residuals, DEM error). Each section also launches the
   random-timeseries comparison routine and saves the figures under its
   respective subdirectory.
3. **Section 4.0 + .i subsections:** After the initial height-map preview of the
   selected points, each “Timeseries Comparisons” subsection regenerates the
   detailed plots using the same set of points. If a point becomes invalid for a
   given layer, a valid fallback pixel is picked automatically, logged, and
   highlighted on the layer-specific height map.
4. **Section 5:** Compares the raw velocity, DEM-error-corrected velocity, and
   short-wavelength velocity layers, enforcing common color scales across the
   relevant products.
5. **Section 6:** Range summary which aggregates the reported ranges to make it
    easy to spot unexpected values across products.

Review the executed notebook to verify:

- All figures rendered correctly (no missing files).
- Fallback WKT and maps exist only when necessary.
- Statistics (RMSE, MAD, R²) fall within expected ranges for each section.

This step helps you spot discrepancies between processing versions quickly and
provides a consolidated artifact (the executed notebook + figures) to archive
with the frame’s QA package.
<br />
<br />
<br />
<br />

## [OPTIONAL] dolphin Software Workflow

Explore the [dolphin workflow README.md and tools](https://github.com/OPERA-Cal-Val/calval-DISP/tree/main/dolphin_PSTworkflow) used to generated sample DISP products with the [dolphin software](https://github.com/isce-framework/dolphin), which is the core software used by OPERA to generate official DISP products.

<br />
<br />
<br />
<br />

## Legacy Tools

Legacy tools can be found under the [legacy subdirectory](https://github.com/OPERA-Cal-Val/calval-DISP/tree/main/legacy).

<br />
<br />
<br />
<br />

## References

Fattahi, H., Simons, M., and Agram, P. (2017). InSAR time-series estimation of the ionospheric phase delay: An extension of the split range-spectrum technique. IEEE Transactions on Geoscience and Remote Sensing, 55(10), 5984-5996, https://doi.org/10.1109/TGRS.2017.2718566

Jolivet, R., Grandin, R., Lasserre, C., Doin, M. P., and Peltzer, G. (2011). Systematic InSAR tropospheric phase delay corrections from global meteorological reanalysis data. Geophysical Research Letters, 38(17), https://doi.org/10.1029/2011GL048757

Doin, M. P., Lasserre, C., Peltzer, G., Cavalié, O., and Doubre, C. (2009). Corrections of stratified tropospheric delays in SAR interferometry: Validation with global atmospheric models. Journal of Applied Geophysics, 69(1), 35-50, https://doi.org/10.1016/j.jappgeo.2009.03.010 

Agnew, D. C. (2012). SPOTL: Some programs for ocean-tide loading.
<br />
<br />
<br />
<br />

## Contributors

M Grace Bato
<br />
Brett Buzzanga
<br />
Marin Govorcin
<br />
Alexander Handwerger
<br />
Jinwoo Kim
<br />
Bryan Raimbault
<br />
Simran Sangha
