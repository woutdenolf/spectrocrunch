#!/bin/bash
# 
# Build PyTMM
# 

# ============Initialize environment============
SCRIPT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $SCRIPT_ROOT/funcs.sh
initEnv

# ============Build PyTMM============
cprint "Download PyTMM ..."
mkdir -p PyTMM
cd PyTMM
if [[ $NOTDRY == true ]]; then
    if [[ ! -d db ]]; then
        git clone https://github.com/polyanskiy/refractiveindex.info-database db
    fi
    if [[ ! -d PyTMM ]]; then
        git clone https://github.com/kitchenknif/PyTMM PyTMM
    fi
fi

cprint "Install PyTMM ..."
if [[ $NOTDRY == true ]]; then
    cd PyTMM

    $PYTHONBIN setup.py build -f
    $PYTHONBIN setup.py install
    
    cd ../db
    cp -R database $PYTHON_PACKAGEDIR/visirlib
fi

BUILDSTEP=$(( $BUILDSTEP+1 ))
BUILDSTEPS=$(( $BUILDSTEPS+1 ))
cd $INSTALL_WD


