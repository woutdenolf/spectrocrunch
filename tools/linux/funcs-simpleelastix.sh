#!/bin/bash
# 
# Install simpleelastix on Linux.
# 

SCRIPT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${SCRIPT_ROOT}/funcs.sh
source ${SCRIPT_ROOT}/funcs-cmake.sh
source ${SCRIPT_ROOT}/funcs-swig.sh
source ${SCRIPT_ROOT}/funcs-python.sh

function simpleelastix_build_dependencies()
{
    local tmp=$(pwd)
    cd ${1}

    mapt-get install lua5.1 liblua5.1-dev
    pip_install virtualenv
    require_cmake 3
    require_swig 3

    cd ${tmp}
}


function simpleelastix_install_fromsource()
{
    local restorewd=$(pwd)

    cprint "Download SimpleElastix ..."
    mkdir -p simpleelastix
    cd simpleelastix
    if [[ $(dryrun) == false && ! -d SimpleElastix ]]; then
        require_web_access
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
            simpleelastix_build_dependencies ${restorewd}

            # http://simpleelastix.readthedocs.io/GettingStarted.html#manually-building-on-linux
            CMAKE_PARAMS="-DBUILD_EXAMPLES:BOOL=OFF \
                          -DBUILD_SHARED_LIBS:BOOL=OFF \
                          -DBUILD_TESTING:BOOL=OFF \
                          -DSimpleITK_USE_SYSTEM_SWIG:BOOL=ON \
                          -DSimpleITK_USE_SYSTEM_LUA:BOOL=ON \
                          -DSimpleITK_USE_SYSTEM_VIRTUALENV:BOOL=ON \
                          -DSimpleITK_USE_SYSTEM_ELASTIX:BOOL=OFF \
                          -DSimpleITK_USE_SYSTEM_ITK:BOOL=OFF \
                          -DPYTHON_EXECUTABLE:FILEPATH=$(python_bin) \
                          -DPYTHON_INCLUDE_DIR:PATH=$(python_include) \
                          -DPYTHON_LIBRARY:FILEPATH=$(python_lib) \
                          -DWRAP_DEFAULT:BOOL=OFF \
                          -DWRAP_PYTHON:BOOL=ON \
                          -DSWIG_DIR:PATH=$(cmd_path swig) \
                          -DSWIG_EXECUTABLE:FILEPATH=$(cmd_full_bin swig)"

            mexec mkdir -p ${prefix}
            cmake -LAH -DCMAKE_INSTALL_PREFIX:PATH="${prefix}" $CMAKE_PARAMS ../SimpleElastix/SuperBuild
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



