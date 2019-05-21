#!/bin/bash
# 
# Install GNU libraries.
# 

SCRIPT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${SCRIPT_ROOT}/funcs-make.sh


function gnu_url()
{
    echo "http://ftp.gnu.org/gnu/${1}/"
}


function gnu_all_versions()
{
    curl -s $(gnu_url ${1}) | grep -E -o ">${1}-[0-9\.]+.tar.gz<" | grep -E -o "[0-9\.]+[0-9]" | sort --version-sort
}


function gnu_latest()
{
    local libname=${1}
    local rversion=${2}
    
    function _gnu_all_versions()
    {
        gnu_all_versions ${libname}
    }
    
    latest_version _gnu_all_versions ${rversion}
}


function gnu_download()
{
    if [[ ! -f ${2}.tar.gz ]]; then
        curl -kL $(gnu_url)/${1}/${2}.tar.gz --output ${2}.tar.gz
    fi
}


function gnu_source_install()
{
    local restorewd=$(pwd)
    local program=${1}
    local rversion=${2}

    cprint "Download ${program} ..."
    mkdir -p ${program}
    cd ${program}

    local version=$(get_local_version)
    if [[ -z ${version} ]]; then
        require_web_essentials
        version=$(gnu_latest ${program} ${rversion})
    fi
    
    local base=${program}-${version}
    echo ${base}
    if [[ $(dryrun) == false && ! -d ${base} ]]; then
        gnu_download ${program} ${base}
    fi

    if [[ $(dryrun) == false ]]; then
        cmake_build_dependencies
        easymake ${program} \
                 ${version}
    fi

    cd ${restorewd}
}


