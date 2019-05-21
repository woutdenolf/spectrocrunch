#!/bin/bash
# 
# Install swig.
# 

SCRIPT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${SCRIPT_ROOT}/funcs-make.sh


function swig_url()
{
    echo "https://sourceforge.net/projects/swig/files/swig/"
}


function swig_all_versions()
{
    curl -ksL $(swig_url) | grep -E -o "title=\"swig-[0-9\.]+" | grep -E -o "[0-9\.]+[0-9]" | sort --version-sort
}


function swig_latest()
{
    latest_version swig_all_versions ${1}
}


function swig_download()
{
    if [[ ! -f ${1}.tar.gz ]]; then
        curl -kL $(swig_url)/${1}/${1}.tar.gz --output ${1}.tar.gz
    fi
}


function swig_build_dependencies()
{
    require_build_essentials
    require_web_essentials
    mapt-get install libpcre3 libpcre3-dev
}


function swig_run_dependencies()
{
    :
}


function swig_source_install()
{
    source_install swig ${1}
}


function swig_version()
{
    if [[ $(swig_exists) == false ]]; then
        echo 0
    else
        swig -version | head -2 | tail -1 | awk '{print $3}'
    fi
}


function swig_exists()
{
    cmdexists swig
}


function swig_system_install()
{
    hash -d swig # swig2/swig3 conflicts
    # Try system installation
    if [[ $(swig_exists) == false ]]; then
        mapt-get install swig
    fi
}


function require_swig()
{
    require_software swig ${1}
}


