After installing dolphin, update the environment to reflect 
mintpy requirements as so:
conda env update --name dolphin-env  --file extra_dependencies.yml

Using the big island of Hawaii as a case example, 
run the dolphin workflow and generate a velocity fit to the 
derived time-series in 3 steps:

1. Process IFGs spanning from 20221001 to 20230701 across a descending track, 
download and use the `esa_world_cover_2021` water mater, download the 
`glo_30` DEM, using the recommended faster `phass` unwrapping option,
and perform some basic parallelization with 2 parallel jobs, 
across 2 parallel tiles and using 4 threads per worker:
/path/to/calval-DISP/dolphin_PSTworkflow/pst_dolphin_workflow.py 
-s 20221001 -e 20230701 -ao '19.0 20.1962 -155.88 -155.014' 
-op 'des' --threadsperworker 4 --nparalleljobs 2 
--nparalleltiles 2 --ntiles 2
--water-mask-file esa_world_cover_2021  --dem-file glo_30
-o pst_output --unwmethod phass

2. Prep prerequisite MintPy inputs needed to generate a velocity field:
/path/to/calval-DISP/dolphin_PSTworkflow/prep_mintpy.py
-m pst_output/static_CSLCs/
-c "pst_output/dolphin_output/stitched_interferograms/*.zeroed.cor.tif"
-u "pst_output/dolphin_output/stitched_interferograms/*.unw.tif"
--geom-dir pst_output/dolphin_output/stitched_interferograms/geometry
--single-reference
--water-mask-file pst_output/dolphin_output/stitched_interferograms/geometry/esa_world_cover_2021_mask.tif
--dem-file pst_output/dolphin_output/stitched_interferograms/geometry/glo_30_DEM.tif
-o mintpy_output

3. Generate velocity fit:
timeseries2velocity.py mintpy_output/timeseries.h5 -o mintpy_output/velocity.h5
