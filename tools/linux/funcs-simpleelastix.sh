#!/bin/bash
# 
# Install simpleelastix on Linux.
# 

SCRIPT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${SCRIPT_ROOT}/funcs.sh
source ${SCRIPT_ROOT}/funcs-cmake.sh
source ${SCRIPT_ROOT}/funcs-swig.sh
source ${SCRIPT_ROOT}/funcs-python.sh
source ${SCRIPT_ROOT}/funcs-lua.sh


function simpleelastix_build_dependencies()
{
    local tmp=$(pwd)
    cd ${1}

    mapt-get install lua5.1 liblua5.1-dev
    #pip_install virtualenv>=13.0 # does not work on rnice
    require_cmake 3.10
    require_swig 3

    cd ${tmp}
}


function simpleelastix_install_fromsource()
{
    if [[ ! -d simpleelastix && ${ARG_SKIPLONG} == true ]]; then
        cprint "Skipping simpleelastix installation"
        return
    fi

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

            local _SYSTEM_SWIG="-DSimpleITK_USE_SYSTEM_SWIG:BOOL=OFF"
            if [[ $(require_new_version $(swig_version) 3) == false ]]; then
                _SYSTEM_SWIG="-DSimpleITK_USE_SYSTEM_SWIG:BOOL=ON \
                              -DSWIG_EXECUTABLE:FILEPATH=$(cmd_full_bin swig) \
                              -DSWIG_DIR:PATH=$(cmd_path swig) "
            fi
            
            local _SYSTEM_LUA="-DSimpleITK_USE_SYSTEM_LUA:BOOL=OFF"
            if [[ $(require_new_version $(lua_version) 5.1) == false ]]; then
                _SYSTEM_LUA="-DSimpleITK_USE_SYSTEM_LUA:BOOL=ON"
            fi
            
            # http://simpleelastix.readthedocs.io/GettingStarted.html#manually-building-on-linux
            CMAKE_PARAMS="-DBUILD_EXAMPLES:BOOL=OFF \
                          -DBUILD_SHARED_LIBS:BOOL=OFF \
                          -DBUILD_TESTING:BOOL=OFF \
                          ${_SYSTEM_SWIG} \
                          ${_SYSTEM_LUA} \
                          -DSimpleITK_USE_SYSTEM_VIRTUALENV:BOOL=OFF \
                          -DSimpleITK_USE_SYSTEM_ELASTIX:BOOL=OFF \
                          -DSimpleITK_USE_SYSTEM_ITK:BOOL=OFF \
                          -DPYTHON_EXECUTABLE:FILEPATH=$(python_bin) \
                          -DPYTHON_INCLUDE_DIR:PATH=$(python_include) \
                          -DPYTHON_LIBRARY:FILEPATH=$(python_lib) \
                          -DWRAP_DEFAULT:BOOL=OFF \
                          -DWRAP_PYTHON:BOOL=ON"

            mexec mkdir -p ${prefix}
            cmake -LAH -DCMAKE_INSTALL_PREFIX:PATH="${prefix}" $CMAKE_PARAMS ../SimpleElastix/SuperBuild
            cmake -DCMAKE_INSTALL_PREFIX:PATH="${prefix}" $CMAKE_PARAMS ../SimpleElastix/SuperBuild
        fi

        cprint "Build SimpleElastix ..."
        if [[ $(dryrun) == false ]]; then
            #OMP_NUM_THREADS=2
            #make -s -j2
            mmakeinstall "simpleelastix"
            #mmakepack "simpleelastix"
            #dpkg -i s ...
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



