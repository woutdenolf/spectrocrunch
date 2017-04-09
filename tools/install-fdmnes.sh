#!/bin/bash
# 
# Install fdmnes on Linux.
# 

# ============Initialize environment============
SCRIPT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $SCRIPT_ROOT/funcs.sh
initEnv

# ============Install fdmnes============
mkdir -p fdmnes
cd fdmnes

if [ ! -d /sware/exp/fdmnes ]; then
    FDMNESVERSION=`curl -s https://api.github.com/repos/gudasergey/fdmnes/releases | grep tag_name | head -n 1 | cut -d '"' -f 4`
    FDMNESZIPNAME=fdmnes_${FDMNESVERSION}.zip
    FDMNESLINK=http://neel.cnrs.fr/IMG/zip/${FDMNESZIPNAME}

    if [ ! -f ${FDMNESZIPNAME} ]; then
        cprint "Download fdmnes ..."

        if [[ $NOTDRY == true ]]; then
            curl -O ${FDMNESLINK}
        fi
    fi
fi

cprint "Install fdmnes ..."
if [[ $NOTDRY == true ]]; then
    # Install in opt
    if [ -d /sware/exp/fdmnes ]; then
        # Link to ESRF sware
        if [ ! -d $SPECTROCRUNCHOPT/fdmnes ]; then
            mexec "ln -s /sware/exp/fdmnes $SPECTROCRUNCHOPT/fdmnes"
        fi
    else
        mexec "unzip -o ${FDMNESZIPNAME} -d $SPECTROCRUNCHOPT"
    fi

    # Link in bin
    if [[ ! -f $SPECTROCRUNCHLOCAL/bin/fdmnes ]]; then
        mexec "mkdir -p $SPECTROCRUNCHLOCAL/bin"
        mexec "ln -s $SPECTROCRUNCHOPT/fdmnes/fdmnes_linux64 $SPECTROCRUNCHLOCAL/bin/fdmnes"
    fi

    # Environment
    addProfile $SPECTROCRUNCHRC "# Installed fdmnes: $SPECTROCRUNCHOPTSTR/fdmnes"
    if [ ! -d /sware/exp/fdmnes ]; then
        addBinPath $SPECTROCRUNCHLOCAL/bin
        addBinPathProfile $SPECTROCRUNCHRC "$SPECTROCRUNCHLOCALSTR/bin"
    fi
fi

cprint "Download pyfdmnes ..."
if [[ $NOTDRY == true && ! -d pyFDMNES ]]; then
    git clone https://github.com/woutdenolf/pyFDMNES.git pyFDMNES
fi

cprint "Install pyfdmnes ..."
if [[ $NOTDRY == true ]]; then
    cd pyFDMNES

    echo "[global]" >> setup.cfg
    echo "fdmnes_path=$SPECTROCRUNCHLOCAL/bin/fdmnes" >> setup.cfg

    $PYTHONBIN setup.py build -f
    $PYTHONBIN setup.py install
fi

BUILDSTEP=$(( $BUILDSTEP+1 ))
BUILDSTEPS=$(( $BUILDSTEPS+1 ))
cd $INSTALL_WD

