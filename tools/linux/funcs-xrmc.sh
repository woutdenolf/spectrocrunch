#!/bin/bash
# 
# Install xrmc.
# 

SCRIPT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${SCRIPT_ROOT}/funcs-make.sh
source ${SCRIPT_ROOT}/funcs-xmimsim.sh


function xrmc_build_dependencies()
{
    require_build_essentials
    require_xmimsim
}


function xrmc_run_dependencies()
{
    require_xmimsim
}


function xrmc_all_versions()
{
    versions_from_github "golosio" "xrmc" "XRMC-[0-9\.]+"
}


function xrmc_latest()
{
    latest_version xrmc_all_versions ${1}
}


function xrmc_download()
{
    if [[ ! -f ${1}.tar.gz ]]; then
        local _rname=${1/xrmc/XRMC}
        curl -L https://github.com/golosio/xrmc/archive/${_rname}.tar.gz --output ${1}.tar.gz
    fi
}


function xrmc_source_install()
{
    #if [[ ! -d xrmc && ${ARG_SKIPLONG} == true ]]; then
    #    cprint "Skipping xrmc installation"
    #    return
    #fi
    
    source_install xrmc "${1}"
}


function xrmc_system_install()
{
    :
}


function xrmc_exists()
{
    cmdexists xrmc
}


function xrmc_version()
{
    if [[ $(xrmc_exists) == false ]]; then
        echo 0
    else
        type xrmc | grep -o -E "/[0-9\.]+[0-9]/" | grep -o -E "[0-9\.]+[0-9]"
    fi
}


function require_xrmc()
{
    require_software xrmc ${1}
}

