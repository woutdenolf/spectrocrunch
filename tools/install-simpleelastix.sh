#!/bin/bash
# 
# Install simpleelastix on Linux.
# 

# ============Initialize environment============
SCRIPT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $SCRIPT_ROOT/funcs.sh
initEnv

# ============Install dependencies============
cprint "Install SimpleElastix dependencies ..."
if [[ $NOTDRY == true && $SYSTEM_PRIVILIGES == true ]]; then
    mexec "apt-get -y install swig"
    mexec "apt-get -y install $PYTHONBINAPT-dev"
fi

# cmake
dpkg --compare-versions "$(cmake --version | head -1 | awk '{print $3}')" "lt" "3.0.0"
if [ $? = 0 ]; then
    source $SCRIPT_ROOT/install-cmake.sh
    cd $INSTALL_WD
fi

# swig
dpkg --compare-versions "$(swig -version | head -2 | tail -1 | awk '{print $3}')" "lt" "3.0.0"
if [ $? = 0 ]; then
    source $SCRIPT_ROOT/install-swig.sh
    cd $INSTALL_WD
fi

# ITK
#source $SCRIPT_ROOT/install-itk.sh
#cd $INSTALL_WD

# ============Install simpleelastix============
if [ ! -f simpleelastix/build/SimpleITK-build/Wrapping/Python/Packaging/setup.py ]; then
    cprint "Download SimpleElastix ..."
    mkdir -p simpleelastix
    cd simpleelastix
    if [[ $NOTDRY == true && ! -d SimpleElastix ]]; then
        git clone https://github.com/kaspermarstal/SimpleElastix SimpleElastix
    fi
    mkdir -p build
    cd build

    if [[ $TIMELEFT == true && ! -f Makefile ]]; then
        cprint "Configure SimpleElastix ..."
        if [[ $NOTDRY == true ]]; then
            # http://simpleelastix.readthedocs.io/GettingStarted.html#manually-building-on-linux
            #
            # TODO:  
            CMAKE_PARAMS="-DBUILD_EXAMPLES:BOOL=OFF \
                          -DBUILD_SHARED_LIBS:BOOL=OFF \
                          -DBUILD_TESTING:BOOL=OFF \
                          -DSIMPLEITK_USE_SYSTEM_SWIG:BOOL=ON \
                          -DSIMPLEITK_USE_SYSTEM_LUA:BOOL=OFF \
                          -DUSE_SYSTEM_VIRTUALENV:BOOL=OFF \
                          -DUSE_SYSTEM_ELASTIX:BOOL=OFF \
                          -DUSE_SYSTEM_ITK:BOOL=OFF \
                          -DPYTHON_EXECUTABLE:FILEPATH=$PYTHON_EXECUTABLE \
                          -DPYTHON_INCLUDE_DIR:PATH=$PYTHON_INCLUDE_DIR \
                          -DPYTHON_LIBRARY:FILEPATH=$PYTHON_LIBRARY \
                          -DITK_DIR:PATH=${ITK_DIR} \
                          -DWRAP_CSHARP:BOOL=OFF \
                          -DWRAP_JAVA:BOOL=OFF \
                          -DWRAP_LUA:BOOL=OFF \
                          -DWRAP_PYTHON:BOOL=ON \
                          -DWRAP_R:BOOL=OFF \
                          -DWRAP_RUBY:BOOL=OFF \
                          -DWRAP_TCL:BOOL=OFF"

            if [[ $INSTALL_SYSTEMWIDE == true ]]; then
                cmake $CMAKE_PARAMS ../SimpleElastix/SuperBuild
            else
                mkdir -p $SPECTROCRUNCHLOCAL
                cmake -DCMAKE_INSTALL_PREFIX:PATH="$SPECTROCRUNCHLOCAL" $CMAKE_PARAMS ../SimpleElastix/SuperBuild
            fi
        fi
        
        BUILDSTEP=$(( $BUILDSTEP+1 ))
    fi

    if [[ $TIMELEFT == true ]]; then
        cprint "Build SimpleElastix ..."
        OMP_NUM_THREADS=2
        if [[ $NOTDRY == true ]]; then
            make -s -j2
            if [[ $TIMELIMITED == true ]]; then
                TIMELEFT=false
            fi
        fi
        BUILDSTEP=$(( $BUILDSTEP+1 ))
    fi
else
    cd simpleelastix/build
    BUILDSTEP=$(( $BUILDSTEP+2 ))
fi

BUILDSTEPS=$(( $BUILDSTEPS+2 ))

if [[ $TIMELEFT == true ]]; then
    cprint "Install SimpleElastix ..."
    if [[ $NOTDRY == true ]]; then
        cd ./SimpleITK-build/Wrapping/Python/Packaging
        $PYTHONBIN setup.py install

        addProfile $SPECTROCRUNCHRC "# Installed simpleelastix: $SPECTROCRUNCHLOCALSTR"
        addLibPath $SPECTROCRUNCHLOCAL/lib
        addLibPathProfile $SPECTROCRUNCHRC "$SPECTROCRUNCHLOCALSTR/lib"
    fi

    BUILDSTEP=$(( $BUILDSTEP+1 ))
fi

BUILDSTEPS=$(( $BUILDSTEPS+1 ))
cd $INSTALL_WD

