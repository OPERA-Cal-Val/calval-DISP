* Note, copy/replace '/u/data-drive/username/' with your own preferred installation path.

Add the following to your `.bashrc`, with a valid path `/u/data-drive/username/scratch` specified.
This is where temporary outputs from dolphin will be written to and managed
export TMPDIR="/u/data-drive/username/scratch"
export TMP=$TMPDIR
export TEMP=$TMPDIR

Rerun bash to refresh your environment

We recommend using the Miniforge conda environment manager, which uses conda-forge as its default code repo
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -p miniforge
miniforge/bin/mamba init bash
conda activate base

First install the main repository dependencies
cd miniforge
mamba env create -f environment.yml
conda activate calval_disp

After initializing the, update the environment to reflect 
other workflow prerequisites as so:

First clone this repo and add it to your path
cd /u/data-drive/username/
git clone https://github.com/OPERA-Cal-Val/calval-DISP.git
setenv PATH ${PATH}:/u/data-drive/username/calval-DISP/dolphin_PSTworkflow

# then install and initialize MintPy build
git clone https://github.com/insarlab/MintPy.git
mamba install --channel conda-forge --file MintPy/requirements.txt
export MINTPY_HOME=/u/data-drive/username/MintPy
export PYTHONPATH=/u/data-drive/username/src
export PATH="${PATH}:${MINTPY_HOME}/src/mintpy/cli"
# then install dolphin-related dependencies
mamba env update --file dolphin_PSTworkflow/extra_dependencies.yml

Initialize `git-lfs` module and re-pull from repo in order to access large CSV file mapping PST S3 database
git lfs install
cd ../
git lfs pull

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
-o mintpy_output
