# OPERA Surface Displacement (DISP) validation tools
Repository of tools intended to address the Cal/Val validation plan for the Level 3 OPERA DISP product suite. Namely, verify that products meet the established OPERA L3 DISP uncertainty requirements.

## Contents

1. [Installation](#installation)
2. [Overview of notebook](#overview-of-notebook)
3. [General Procedures](#general-procedures)
    -   [InSAR LOS velocity estimation](#insar-los-velocity-estimation)
    -   [Validation data LOS velocity estimation](#validation-data-los-velocity-estimation)
    -   [Statistical analyses](#statistical-analyses)
4. [Validation Approaches](#validation-approaches)
    -   [Validation approach 1 (VA1)](#validation-approach-1)
    -   [Validation approach 2 (VA2)](#validation-approach-2)
    -   [Validation approach 3 (VA3)](#validation-approach-3)
5. [References](#references)
6. [Contributors](#contributors)

## Installation

A detailed installation on how to build a working environment can be found __[here](https://github.com/OPERA-Cal-Val/calval-DISP/blob/main/docs/installation.md/)__.

## Overview of notebook

DISP analog, standardized interferometric products can be obtained from the Alaska Satellite Facility (ASF) Hybrid Pluggable Processing Pipeline (HyP3) on demand service, and in the form of the Advanced Rapid Image Analysis (ARIA) project’s formulated Geocoded Unwrapped Phase product (ARIA-GUNW).

Supporting input of such stacks of analog products across a range of diverse study areas, the __[DISP-S1_Requirement_Validation.ipynb](https://github.com/OPERA-Cal-Val/calval-DISP/blob/main/DISP-S1_Requirement_Validation.ipynb/)__ notebook contains a suite of functionality to evaluate the uncertainties associated with the validation requirements as outlined in the OPERA validation plan. The notebook itself was adopted the NISAR team Algorithm Theoretical Basis Document (ATBD) __[notebook for secular motion](https://github.com/nisar-solid/ATBD/blob/main/methods/secular/Secular_Requirement_Validation.ipynb/)__.

## General Procedures

Summarized below are the general procedures to implement our validation approaches (outlined in the next section) including the two primary statistical analyses that are performed:
1.	InSAR LOS velocity estimation: We estimate the long term line-of-sight rates and corresponding uncertainties by fitting a time-series function which includes a linear rate parameter across the displacement history, derived from the reconstructed displacement time-series using consecutive DISP products in time. Where needed, we may include additional complexities (e.g. Heaviside step, exponential or logarithmic function for coseismic and postseismic events/processes).
2.	Validation data LOS velocity estimation: We project the independent Validation displacement data (e.g GNSS time-series) in the radar look direction and then estimate the velocity over the same period of time and with the same time-series models as from step 1. The GNSS time-series data will have its processing and noise corrections applied. As an example, the GNSS time-series data obtained from the University of Nevada, Reno (UNR) is corrected for instrument changes, ionospheric and tropospheric propagation delays, as well as solid earth tides [Blewitt et al., 2018].
3.	Statistical analyses:  We mainly implement similar analyses as exemplified to assess if the OPERA L3 DISP product requirements are achieved.
  - a)	Variogram analysis for the line-of-sight velocity residuals will be performed. Using the structure function, values at length scales 0.1 km < L < 50 km will be extracted and then compared to the OPERA L3 DISP requirements for validation.
  - b)	We calculate the mean and the standard deviation of the velocity residuals (e.g. perform scatter plot/RMSE analysis) and then assess if the OPERA L3 DISP requirements are satisfied.

For all the validation analyses, we report the accuracy that can be achieved with and without applying the optional correction layers embedded within the L3 DISP products, including the ionospheric propagation delays (e.g., [Fattahi et al., 2017]), tropospheric propagation delays (e.g., [Jolivet et al., 2011] and [Doin et al., 2009]), and solid earth tide effects (e.g. [Agnew 2012]).
<br />
<br />
<br />
<br />


## Validation Approaches

Three validation approaches will be performed separately depending on the Validation dataset that will be used to assess the DISP product requirements. Each validation approach is dubbed in this document as Validation Approach n (VAn) for n = 1, 2, 3, i.e. VA1, VA2, and VA3. The validation approaches are primarily differentiated based on the surface displacement occurring at the validation site (i.e. deforming or negligibly deforming) and/or the instrument used to validate the InSAR measurement (i.e. GNSS, corner reflectors, InSAR). NOTE: also supported are select analyses involving the Fine Resolution InSAR using Generalized Eigenvectors (FRInGE) __[package](https://github.com/isce-framework/fringe/)__, which directly estimates deformation time-series by exploiting interferometric covariance matrices such that intermediate interferometric products are not generated.

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


## References

Fattahi, H., Simons, M., and Agram, P. (2017). InSAR time-series estimation of the ionospheric phase delay: An extension of the split range-spectrum technique. IEEE Transactions on Geoscience and Remote Sensing, 55(10), 5984-5996, https://doi.org/10.1109/TGRS.2017.2718566

Jolivet, R., Grandin, R., Lasserre, C., Doin, M. P., and Peltzer, G. (2011). Systematic InSAR tropospheric phase delay corrections from global meteorological reanalysis data. Geophysical Research Letters, 38(17), https://doi.org/10.1029/2011GL048757

Doin, M. P., Lasserre, C., Peltzer, G., Cavalié, O., and Doubre, C. (2009). Corrections of stratified tropospheric delays in SAR interferometry: Validation with global atmospheric models. Journal of Applied Geophysics, 69(1), 35-50, https://doi.org/10.1016/j.jappgeo.2009.03.010 

Agnew, D. C. (2012). SPOTL: Some programs for ocean-tide loading.

## Contributors

M Grace Bato
<br />
Brett Buzzanga
<br />
Marin Govorcin
<br />
Alexander Handwerger
<br />
Simran Sangha
