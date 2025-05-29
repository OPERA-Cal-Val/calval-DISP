## Install OPERA DISP Cal/Val environment

Instructions derived and modified from https://github.com/nisar-solid/ATBD/blob/main/docs/installation.md 

## 1. Install conda

### Initial environment setup
* Note, copy/replace '/u/data-drive/username/' with your own preferred installation path.

### Setup conda
We recommend using the Miniforge conda environment manager, which uses conda-forge as its default code repo:
```bash
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -p miniforge
miniforge/bin/mamba init bash
conda activate base
```

Rerun bash to refresh your environment:
```bash
bash
```

### Install main repository dependencies
```bash
cd /u/data-drive/username/
git clone https://github.com/OPERA-Cal-Val/calval-DISP.git
mamba env create -f calval-DISP/environment.yml
conda activate calval_disp

# add repo tools to your path
export PATH="${PATH}:/u/data-drive/username/calval-DISP"

# set cap to circumvent potential dolphin crash
export XLA_PYTHON_CLIENT_MEM_FRACTION=".10"
```

## 2. Install other prerequisite packages

### Install development version of MintPy
```bash
git clone https://github.com/insarlab/MintPy.git
mamba update --file MintPy/requirements.txt --name calval_disp
export MINTPY_HOME=/u/data-drive/username/MintPy
export PYTHONPATH=/u/data-drive/username/src
export PATH="${PATH}:${MINTPY_HOME}/src/mintpy/cli"
```

### Install development version of RAiDER
```bash
git clone https://github.com/dbekaert/RAiDER.git
cd RAiDER
mamba env update --file environment.yml --name calval_disp
python -m pip install -e .
cd ../
```

### Set PYTHONPATH
```bash
# set paths to prerequisite tools
export PYTHONPATH="${PYTHONPATH}:/u/data-drive/username/calval-DISP"
cd ../
```

### Install papermill and initiate kernel

```bash
python3 -m pip install papermill
python -m ipykernel install --user --name calval_disp
```

### Test the installation

Run the following to test the installation:

```bash
# Display help for the download script (try using 'python' if issues occur)
run1_download_DISP_S1_Static.py --h 

# Display help for MintPy
smallbaselineApp.py -h
```