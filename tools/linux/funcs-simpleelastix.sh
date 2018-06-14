#!/bin/bash
# 
# Install simpleelastix on Linux.
# 

SCRIPT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $SCRIPT_ROOT/funcs.sh
source $SCRIPT_ROOT/funcs-cmake.sh
source $SCRIPT_ROOT/funcs-swig.sh

function simpleelastix_build_dependencies()
{
    require_cmake 3
    require_swig 3
}


function simpleelastix_install_fromsource()
{
    local restorewd=$(pwd)

    cprint "Download SimpleElastix ..."
    mkdir -p simpleelastix
    cd simpleelastix
    if [[ $(dryrun) == false && ! -d SimpleElastix ]]; then
        requires_web_access
        git clone https://github.com/kaspermarstal/SimpleElastix SimpleElastix
    fi
    mkdir -p build
    cd build

    local prefix=$(project_opt)/simpleelastics
    local prefixstr=$(project_optstr)/simpleelastics
    local outdir=./SimpleITK-build/Wrapping/Python/Packaging
    if [[ $(dryrun) == false && ! -f ${outdir}/setup.py ]]; then

        cprint "Configure SimpleElastix ..."
        if [[ $(dryrun) == false ]]; then
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
                          -DPYTHON_EXECUTABLE:FILEPATH=$(python_bin) \
                          -DPYTHON_INCLUDE_DIR:PATH=$(python_include) \
                          -DPYTHON_LIBRARY:FILEPATH=$(python_lib) \
                          -DWRAP_CSHARP:BOOL=OFF \
                          -DWRAP_JAVA:BOOL=OFF \
                          -DWRAP_LUA:BOOL=OFF \
                          -DWRAP_PYTHON:BOOL=ON \
                          -DWRAP_R:BOOL=OFF \
                          -DWRAP_RUBY:BOOL=OFF \
                          -DWRAP_TCL:BOOL=OFF"
                          #-DITK_DIR:PATH=${ITK_DIR}

            simpleelastix_build_dependencies
            mkdir -p ${prefix}
            cmake -DCMAKE_INSTALL_PREFIX:PATH="${prefix}" $CMAKE_PARAMS ../SimpleElastix/SuperBuild
        fi
        
        cprint "Build SimpleElastix ..."
        if [[ $(dryrun) == false ]]; then
            OMP_NUM_THREADS=2
            make -s -j2
        fi
    fi

    cprint "Install SimpleElastix ..."
    if [[ $(dryrun) == false ]]; then
        cd ${outdir}
        $(python_bin) setup.py install

        addProfile $(project_resource) "# Installed simpleelastix: ${prefixstr}"
        addLibPath ${prefix}/lib
        addLibPathProfile $(project_resource) "${prefixstr}/lib"
    fi

    cd ${restorewd}
}



function require_simpleelastix()
{
    cprintstart
    cprint "Verify simpleelastix ..."

    # Requirements (for running)
    require_python

    # Check
    if [[ $(python_hasmodule "SimpleITK") == true ]]; then
        cprint "Python module \"simpleelastix\" is installed"
        cprintend
        return
    fi

    # Install from source
    simpleelastix_install_fromsource

    # Check
    if [[ $(python_hasmodule "SimpleITK") == true ]]; then
        cprint "Python module \"simpleelastix\" is installed"
    else
        cprint "Python module \"simpleelastix\" is NOT installed"
    fi

    cprintend
}



