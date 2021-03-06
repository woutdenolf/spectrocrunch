#!/bin/bash
# 
# Install libgeos on Linux.
# 

# ============Initialize environment============
SCRIPT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $SCRIPT_ROOT/funcs.sh
initEnv

# ============Install dependencies============
cprint "Install libgeos dependencies ..."

# cmake
dpkg --compare-versions "$(cmake --version | head -1 | awk '{print $3}')" "lt" "3.1.3"
if [ $? = 0 ]; then
    source $SCRIPT_ROOT/install-cmake.sh
    cd $INSTALL_WD
fi

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
        PREFIX=$(project_prefix)
        PREFIXSTR=$(project_prefix)
        mkdir -p ${PREFIX}
        cmake -DCMAKE_INSTALL_PREFIX:PATH="${PREFIX}" ../geos
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
    mmakeinstall

    addProfile $(project_resource) "# Installed libgeos: ${PREFIXSTR}"
    addLibPath ${PREFIX}/lib
    addLibPathProfile $(project_resource) "${PREFIXSTR}/lib"
fi

BUILDSTEP=$(( $BUILDSTEP+1 ))
BUILDSTEPS=$(( $BUILDSTEPS+1 ))
cd $INSTALL_WD

