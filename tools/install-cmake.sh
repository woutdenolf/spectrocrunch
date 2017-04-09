#!/bin/bash
# 
# Install xraylib on Linux.
# 

# ============Initialize environment============
SCRIPT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $SCRIPT_ROOT/funcs.sh
initEnv

# ============Install cmake============
cprint "Download cmake ..."
mkdir -p cmake
cd cmake
if [[ $NOTDRY == true && ! -d cmake-3.7.2 ]]; then
    wget --no-check-certificate http://www.cmake.org/files/v3.7/cmake-3.7.2.tar.gz
    tar -xzf cmake-3.7.2.tar.gz
fi
cd cmake-3.7.2

if [[ ! -f Makefile ]]; then
    cprint "Configure cmake ..."
    if [[ $NOTDRY == true ]]; then
        if [[ $SYSTEM_PRIVILIGES == true ]]; then
            ./configure
        else
            mkdir -p $SPECTROCRUNCHLOCAL
            ./configure --prefix=$SPECTROCRUNCHLOCAL
        fi
    fi
    cprint "Build cmake ..."
    if [[ $NOTDRY == true ]]; then
        make -s -j2
    fi
fi


cprint "Install cmake ..."
if [[ $NOTDRY == true ]]; then
    mexec "make install -s"

    # Add path just for this installation script
    addBinPath $SPECTROCRUNCHLOCAL/bin
fi

BUILDSTEP=$(( $BUILDSTEP+1 ))
BUILDSTEPS=$(( $BUILDSTEPS+1 ))
cd $INSTALL_WD

