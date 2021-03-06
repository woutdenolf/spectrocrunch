#!/bin/bash
# 
# Install simpleelastix on Linux.
# 

SCRIPT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${SCRIPT_ROOT}/funcs-make.sh
source ${SCRIPT_ROOT}/funcs-cmake.sh
source ${SCRIPT_ROOT}/funcs-swig.sh
source ${SCRIPT_ROOT}/funcs-python.sh
source ${SCRIPT_ROOT}/funcs-lua.sh


function simpleelastix_build_dependencies()
{
    # Not sure why we need Lua when we don't need the Lua wrapping
    mapt-get install lua5.3 liblua5.3-dev
    require_pythondev
    require_cmake 3.10
    require_swig 3
}


function simpleelastix_download()
{
    git clone https://github.com/kaspermarstal/SimpleElastix ${1}

    cd ${1}
    
    #git reset --hard 2a79d15
    #git checkout v1.2.4  #Nov 13, 2019
    #git checkout v1.1.0
    #git reset --hard 49af818 # Sep 10, 2018  (v1.1.0)
    #git reset --hard cf75ff4 # Dec 21, 2017  (Git ITK v4.13.0)
    #git reset --hard 3de935c # Oct 30, 2017 
    #git reset --hard v1.1.0
    
    # works from cmake 3.1 onwards:
    #sed -i '/EXCLUDE_FROM_MAIN 1/d' SuperBuild/SuperBuild.cmake

    cd ..
}


function simpleelastix_configure()
{
    if [[ -e "Makefile" ]]; then
        cprint "Configure ${1} (${2}): already configured."
    else
        # Show
        cmake -LAH "${@:3}"
        # Configure and build in one go
        cmake "${@:3}"
    fi
}


function simpleelastix_source_install()
{
    if [[ ! -d simpleelastix && ${ARG_SKIPLONG} == true ]]; then
        cprint "Skipping simpleelastix installation"
        return
    fi

    # Will be called again by "source_install" but
    # we need it already here to modify the arguments:
    simpleelastix_build_dependencies

    local _SYSTEM_SWIG="-DSimpleITK_USE_SYSTEM_SWIG:BOOL=OFF"
    if [[ $(require_new_version $(swig_version) 3) == false ]]; then
        _SYSTEM_SWIG="-DSimpleITK_USE_SYSTEM_SWIG:BOOL=ON \
                      -DSWIG_EXECUTABLE:FILEPATH=$(cmd_full_bin swig) \
                      -DSWIG_DIR:PATH=$(cmd_path swig)"
    fi

    local _SYSTEM_LUA="-DSimpleITK_USE_SYSTEM_LUA:BOOL=OFF"
    if [[ $(require_new_version $(lua_version) 5.3) == false ]]; then
        _SYSTEM_LUA="-DSimpleITK_USE_SYSTEM_LUA:BOOL=ON"
    fi

    local _SYSTEM_VIRTUALENV="-DSimpleITK_USE_SYSTEM_VIRTUALENV:BOOL=OFF"
    if [[ $(python_hasmodule virtualenv) == true || $(python_hasmodule venv) == true ]]; then
        _SYSTEM_VIRTUALENV="-DSimpleITK_USE_SYSTEM_VIRTUALENV:BOOL=ON"
    fi
    
    # http://simpleelastix.readthedocs.io/GettingStarted.html#manually-building-on-linux
    local rversion=${1}
    source_install simpleelastix "${rversion}" \
                  -DCMAKE_INSTALL_PREFIX:PATH='${prefix}' \
                  -DBUILD_EXAMPLES:BOOL=OFF \
                  -DBUILD_SHARED_LIBS:BOOL=OFF \
                  -DBUILD_TESTING:BOOL=OFF \
                  -DBUILD_DOCUMENTS:BOOL=OFF \
                  ${_SYSTEM_SWIG} \
                  ${_SYSTEM_LUA} \
                  -DSimpleITK_PYTHON_USE_VIRTUALENV:BOOL=OFF \
                  -DSimpleITK_USE_SYSTEM_ELASTIX:BOOL=OFF \
                  -DSimpleITK_USE_SYSTEM_ITK:BOOL=OFF \
                  -DSimpleITK_Module_AnisotropicDiffusionLBR:BOOL=ON \
                  -DPYTHON_EXECUTABLE:FILEPATH=$(python_full_bin) \
                  -DPYTHON_INCLUDE_DIR:PATH=$(python_include) \
                  -DPYTHON_LIBRARY:FILEPATH=$(python_lib) \
                  -DWRAP_DEFAULT:BOOL=OFF \
                  -DWRAP_PYTHON:BOOL=ON \
                  -DWRAP_RUBY:BOOL=OFF \
                  -DWRAP_JAVA:BOOL=OFF \
                  -DWRAP_CSHARP:BOOL=OFF \
                  -DWRAP_LUA:BOOL=OFF \
                  -DWRAP_TCL:BOOL=OFF \
                  -DWRAP_R:BOOL=OFF \
                  ../simpleelastix-'${version}'/SuperBuild
}


function simpleelastix_install()
{
    local outdir="SimpleITK-build/Wrapping/Python"
    if [[ -f ${outdir}/Packaging/setup.py ]]; then
        cd ${outdir}
        $(python_bin) Packaging/setup.py install
    else
        cerror "Missing file python ${outdir}/Packaging/setup.py"
    fi
}


function simpleelastix_exists()
{
    python_hasmodule "SimpleITK"
}


function simpleelastix_version()
{
    if [[ $(simpleelastix_exists) == false ]]; then
        echo 0
    else
        python -c "import SimpleITK;print(SimpleITK.Version_VersionString())"
    fi
}


function require_simpleelastix()
{
    require_software simpleelastix
}


