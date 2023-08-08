## Install OPERA DISP Cal/Val environment

Instructions derived and modified from https://github.com/nisar-solid/ATBD/blob/main/docs/installation.md 

### 1. Install conda

```bash
mkdir -p ~/tools; cd ~/tools

# download, install and setup (mini/ana)conda
# for Linux, use Miniconda3-latest-Linux-x86_64.sh
# for macOS, opt 2: curl https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -o Miniconda3-latest-MacOSX-x86_64.sh
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
bash Miniconda3-latest-MacOSX-x86_64.sh -b -p ~/tools/miniconda3
~/tools/miniconda3/bin/conda init bash
```

Close and restart the shell for changes to take effect.

```bash
conda config --add channels conda-forge
conda install git mamba --yes
```

### 2. Install ATBD to `atbd` environment

#### Download source code

```bash
cd ~/tools
git clone https://github.com/nisar-solid/ATBD.git
git clone https://github.com/aria-tools/ARIA-tools.git
```

#### Create `atbd` environment and install pre-requisites

```bash
# create new environment
# install dependencies with mamba by using the local `environment.yml` file
mamba env create -f environment.yml
```

#### Source your installation

Create an alias `load_atbd` in `~/.bash_profile` file for easy activation, _e.g._:

```bash
alias load_atbd='conda activate atbd; source ~/tools/ATBD/docs/config.rc'
```

#### Install MintPy from source (required to access stable version)

Create an alias `load_atbd` in `~/.bash_profile` file for easy activation, _e.g._:

```bash
git clone https://github.com/insarlab/MintPy.git
python -m pip install -e MintPy
```

#### Install papermill and initiate kernel

```bash
python3 -m pip install papermill
python -m ipykernel install --user --name atbd
```

#### Test the installation

Run the following to test the installation:

```bash
load_atbd
ariaDownload.py -h
smallbaselineApp.py -h
```

### 3. Prep Notebook

Open `DISP-S1_Requirement_Validation.ipynb` and add new section in the `sites` dictionary found in the 3rd cell
which conforms to the parameters of your particular AO (e.g. unique name, geographic region, start/end date, reference GNSS station, etc).

Also note that these parameters for particular region of interest may have already been established in the source notebook. So check that our first,
making sure to of course adjust your start/end dates to ensure your series spans 4 years as per the requirement.


### 4. Run Notebook

For example, here is a sample run for the Central Valley, California case study for descending Sentinel-1 track 144.

```bash
#Args:
# the name of your output notebook
# site name (defined within the site dictionary of the input notebook)
# output work directory
papermill DISP-S1_Requirement_Validation.ipynb MYOUTPUT_DISP-S1_Requirement_Validation.ipynb -p site 'CentralValleyD144_4yr' -p work_dir 'D144' -k atbd
```
