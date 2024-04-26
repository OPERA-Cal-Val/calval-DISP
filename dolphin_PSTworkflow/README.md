# OPERA Project Science Team (PST) dolphin workflow 
Repository of tools to leverage the dolphin software to produce pseudo- OPERA Surface Displacement (DISP) interferometric products from OPERA Coregistered Single-Look Complex (CSLC) products used for validation efforts.

------
## Contents

1. [Installation](#installation)
    -   [Initial environment setup](#initial-environment-setup)
    -   [Setup conda](#setup-conda)
    -   [Install main repository dependencies](#install-main-repository-dependencies)
    -   [Install development version of MintPy](#install-development-version-of-MintPy)
    -   [Clone ATBD repo](#clone-ATBD-repo)
2. [Running dolphin workflow](#running-dolphin-workflow)
    -   [Run PST dolphin tool](#run-pst-dolphin-tool)
    -   [Run MintPy preparation tool](#run-MintPy-preparation-tool)
3. [Running validation notebook](#running-validation-notebook)

------
## Installation

### Initial environment setup
* Note, copy/replace '/u/data-drive/username/' with your own preferred installation path.

Add the following to your `.bashrc`, with a valid path `/u/data-drive/username/scratch` specified.
This is where temporary outputs from dolphin will be written to and managed:
```.bash
export TMPDIR="/u/data-drive/username/scratch"
export TMP=$TMPDIR
export TEMP=$TMPDIR
```

Rerun bash to refresh your environment:
```.bash
bash
```

### Setup conda
We recommend using the Miniforge conda environment manager, which uses conda-forge as its default code repo:
```.bash
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -p miniforge
miniforge/bin/mamba init bash
conda activate base
```

### Install main repository dependencies
```.bash
cd /u/data-drive/username/
git clone https://github.com/OPERA-Cal-Val/calval-DISP.git
mamba env create -f calval-DISP/environment.yml
conda activate calval_disp

# install extra dolphin workflow-related dependencies
mamba env update --file calval-DISP/dolphin_PSTworkflow/extra_dependencies.yml --name calval_disp

# add repo tools to your path
export PATH=${PATH}:/u/data-drive/username/calval-DISP/dolphin_PSTworkflow"
```

### Install development version of dolphin
```.bash
###NOTE temp redirect to PR for UTM coordinate bug fix
git clone https://github.com/isce-framework/dolphin.git
cd dolphin
mamba update --name calval_disp --file requirements.txt
python -m pip install .
cd ../
```

### Install development version of MintPy
```.bash
git clone https://github.com/insarlab/MintPy.git
cd MintPy
###NOTE temp redirect to stable revision before major code overhaul
git reset --hard  f3324b8
cd ../
mamba update --name calval_disp --file MintPy/requirements.txt
export MINTPY_HOME=/u/data-drive/username/MintPy
export PYTHONPATH=/u/data-drive/username/src
export PATH="${PATH}:${MINTPY_HOME}/src/mintpy/cli"
```

### Clone ATBD repo
```.bash
###NOTE temp redirect to PR for UTM coordinate bug fix
git clone https://github.com/sssangha/ATBD.git
cd ATBD
git checkout sss_transectnans
# set paths to prerequisite tools
export PYTHONPATH"${PYTHONPATH}:/u/data-drive/username/ATBD"
cd ../
```

------
## Running dolphin workflow

Using the big island of Hawaii as a case example, run the dolphin workflow and generate a velocity fit to the derived time-series in 2 steps:

### Run PST dolphin tool

Process IFGs spanning from 20221001 to 20230701 across a descending track, download and use the `esa_world_cover_2021` water mater, download the `glo_30` DEM, using the recommended faster `phass` unwrapping option, and perform some basic parallelization with 2 parallel jobs, across 2 parallel tiles and using 4 threads per worker:

```.bash
pst_dolphin_workflow.py 
-s 20221001 -e 20230701 -ao '19.0 20.1962 -155.88 -155.014' 
-op 'des' --threadsperworker 4 --nparalleljobs 2 
--nparalleltiles 2 --ntiles 2
--water-mask-file esa_world_cover_2021  --dem-file glo_30
-o pst_output --unwmethod phass
```

### Run MintPy preparation tool

Prep prerequisite MintPy inputs needed to generate a velocity field:
```.bash
prep_mintpy.py
-m pst_output/static_CSLCs/
-c "pst_output/dolphin_output/stitched_interferograms/*.zeroed.cor.tif"
-u "pst_output/dolphin_output/stitched_interferograms/*.unw.tif"
--geom-dir pst_output/dolphin_output/stitched_interferograms/geometry
--single-reference
-o mintpy_output
```

------
## Running validation notebook

You can verify whether your case study meets the DISP validation requirements by running the `DISP-S1_dolphin_Requirement_Validation.ipynb`

First copy over the source notebook found here `/path/to/calval-DISP/dolphin_PSTworkflow/` to your local working directory:
```.bash
cp /u/data-drive/username/calval-DISP/dolphin_PSTworkflow/DISP-S1_dolphin_Requirement_Validation.ipynb .
```

In the 3rd cell, make sure to specify the valid path to your dolphin and mintpy output directories, if they exist. 

Then add an entry or modify and existing entry for a given case study in the 4th `List of CalVal Sites:` which details e.g. the region of interest, time span, reference GNSS station. Note you can rerun the above dolphin workflow from scratch if you specify `use_staged_data` = False.

Finally, run the notebook via papermill, specifying the appropriate `site name` in the dictionary as an argument.
```.bash
papermill DISP-S1_dolphin_Requirement_Validation.ipynb run_D087_DISP-S1_Requirement_Validation.ipynb -p site 'des_D087'
```
