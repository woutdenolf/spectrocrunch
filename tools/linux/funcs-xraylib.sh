#!/bin/bash
# 
# Install xraylib.
# 

SCRIPT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${SCRIPT_ROOT}/funcs-make.sh
source ${SCRIPT_ROOT}/funcs-swig.sh
source ${SCRIPT_ROOT}/funcs-python.sh


function xraylib_build_dependencies()
{
    local tmp=$(pwd)
    cd ${1}

    require_build_essentials
    require_pythondev
    require_swig 3
    pip_install numpy
    pip_install cython

    #$(pip_bin) freeze | grep numpy > requirements.txt
    #pip_install -r requirements.txt
    #rm requirements.txt

    cd ${tmp}
}


function xraylib_run_dependencies()
{
    require_python
}


function xraylib_all_versions()
{
    versions_from_github "tschoonj" "xraylib" "xraylib-[0-9\.]+"
}


function xraylib_latest()
{
    latest_version xraylib_all_versions ${1}
}


function xraylib_download()
{
    if [[ ! -f ${1}.tar.gz ]]; then
        curl -L https://github.com/tschoonj/xraylib/archive/${1}.tar.gz --output ${1}.tar.gz
    fi
}


function xraylib_source_install()
{
    #if [[ ! -d xraylib && ${ARG_SKIPLONG} == true ]]; then
    #    cprint "Skipping xraylib installation"
    #    return
    #fi
    
    source_install xraylib "${1}" \
         --enable-python \
         --enable-python-integration \
         --disable-java \
         --disable-lua \
         --disable-ruby \
         --disable-php \
         --disable-pascal \
         --disable-idl \
         --disable-perl \
         --disable-fortran2003 \
         PYTHON="$(python_full_bin)"
}


function xraylib_system_install()
{
    :
}


function xraylib_exists()
{
    python_hasmodule xraylib
}


function xraylib_version()
{
    if [[ $(xraylib_exists) == false ]]; then
        echo 0
    else
        python -c "import xraylib;print(xraylib.__version__)"
    fi
}


function require_xraylib()
{
    require_software xraylib ${1}
}

