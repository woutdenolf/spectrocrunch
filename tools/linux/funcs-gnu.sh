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


function gnu_download()
{
    if [[ ! -f ${2}.tar.gz ]]; then
        curl -kL $(gnu_url)/${1}/${2}.tar.gz --output ${2}.tar.gz
    fi
}


function glibc_exists()
{
    cmdexists ldd
}


function glibc_version()
{
    if [[ $(glibc_exists) == false ]]; then
        echo 0
    else
        ldd --version | head -1 | grep -E -o '[\.0-9]+[0-9]$'
    fi
}


function require_gnu()
{
    local program=${1}
    local rversion=${2}
    eval "function ${program}_all_versions(){ gnu_all_versions ${program};}"
    eval "function ${program}_latest(){ latest_version ${program}_all_versions \${1};}"
    eval "function ${program}_download(){ gnu_download ${program} \${1};}"
    eval "function ${program}_source_install(){ gnu_source_install ${program} \${1};}"
    require_software ${program} "${rversion}"
}

