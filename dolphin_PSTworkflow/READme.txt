After installing dolphin, update the environment to reflect mintpy requirements as so:
conda env update --name dolphin-env  --file extra_dependencies.yml


Install mintpy as so
cd /PATH/TO/DOLPHIN/CLONE
git clone https://github.com/insarlab/MintPy.git

Add below in your source file, e.g. ~/.bash_profile for bash users or ~/.cshrc for csh/tcsh users:
if [ -z ${PYTHONPATH+x} ]; then export PYTHONPATH=""; fi
export MINTPY_HOME=/PATH/TO/DOLPHIN/CLONE/MintPy
export PATH=${PATH}:${MINTPY_HOME}/src/mintpy/cli
export PYTHONPATH=${PYTHONPATH}:${MINTPY_HOME}/src
