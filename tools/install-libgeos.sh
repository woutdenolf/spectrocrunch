#!/bin/bash
# 
# Install libgeos on Linux.
# 

# ============Initialize environment============
SCRIPT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $SCRIPT_ROOT/funcs.sh
initEnv

# ============Install libgeos============
if [ ! -f libgeos/build/lib/libgeos.so ]; then
    cprint "Download libgeos ..."
    mkdir -p libgeos
    cd libgeos
    if [[ $NOTDRY == true && ! -d geos ]]; then
        git clone https://github.com/OSGeo/geos
    fi
    mkdir -p build
    cd build

    cprint "Configure libgeos ..."
    if [[ $NOTDRY == true ]]; then
        if [[ $SYSTEM_PRIVILIGES == true ]]; then
            cmake ../geos
        else
            mkdir -p $SPECTROCRUNCHLOCAL
            cmake -DCMAKE_INSTALL_PREFIX:PATH="$SPECTROCRUNCHLOCAL" ../geos
        fi
    fi

    cprint "Build libgeos ..."
    if [[ $NOTDRY == true ]]; then
        make -s -j2
        make check
    fi
else
    cd libgeos/build
fi

cprint "Install libgeos ..."
if [[ $NOTDRY == true ]]; then
    mexec "make install -s"

    addProfile $SPECTROCRUNCHRC "# Installed libgeos: $SPECTROCRUNCHLOCALSTR"
    addLibPath $SPECTROCRUNCHLOCAL/lib
    addLibPathProfile $SPECTROCRUNCHRC "$SPECTROCRUNCHLOCALSTR/lib"
fi

BUILDSTEP=$(( $BUILDSTEP+1 ))
BUILDSTEPS=$(( $BUILDSTEPS+1 ))
cd $INSTALL_WD

