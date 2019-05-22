#!/bin/bash
# 
# Install xmimsim.
# 

SCRIPT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${SCRIPT_ROOT}/funcs-make.sh
source ${SCRIPT_ROOT}/funcs-xraylib.sh


function libxml2_url()
{
    echo "ftp://xmlsoft.org/libxml2/"
}


function libxslt_url()
{
    echo "ftp://xmlsoft.org/libxslt/"
}


function libxml2_all_versions()
{
    versions_from_site $(libxml2_url) "libxml2-[0-9\.]+[0-9]\.tar\.gz$"
}


function libxslt_all_versions()
{
    versions_from_site $(libxslt_url) "libxslt-[0-9\.]+[0-9]\.tar\.gz$"
}


function easyrng_all_versions()
{
    versions_from_github "tschoonj" "easyRNG" "easyRNG-[0-9\.]+"
}


function xmimsim_all_versions()
{
    versions_from_github "tschoonj" "xmimsim" "XMI-MSIM-[0-9\.]+"
}


function libxml2_latest()
{
    latest_version libxml2_all_versions ${1}
}


function libxslt_latest()
{
    latest_version libxslt_all_versions ${1}
}


function easyrng_latest()
{
    latest_version easyrng_all_versions ${1}
}


function xmimsim_latest()
{
    latest_version xmimsim_all_versions ${1}
}


function libxml2_download()
{
    if [[ ! -f ${1}.tar.gz ]]; then
        curl -L $(libxml2_url)${1}.tar.gz --output ${1}.tar.gz
    fi
}


function libxslt_download()
{
    if [[ ! -f ${1}.tar.gz ]]; then
        curl -L $(libxslt_url)${1}.tar.gz --output ${1}.tar.gz
    fi
}


function easyrng_download()
{
    if [[ ! -f ${1}.tar.gz ]]; then
        local _rname=${1/easyrng/easyRNG}
        curl -L https://github.com/tschoonj/easyRNG/archive/${_rname}.tar.gz --output ${1}.tar.gz
    fi
}


function xmimsim_download()
{
    if [[ ! -f ${1}.tar.gz ]]; then
        local _rname=${1/xmimsim/XMI-MSIM}
        curl -L https://github.com/tschoonj/xmimsim/archive/${_rname}.tar.gz --output ${1}.tar.gz
    fi
}


function libxml2_build_dependencies()
{
    require_build_essentials
}


function libxslt_build_dependencies()
{
    require_build_essentials
    require_libxml2
}


function easyrng_build_dependencies()
{
    require_build_essentials
}


function xmimsim_build_dependencies()
{
    require_build_essentials
    require_easyrng
    require_libxslt
    require_xraylib
}


function libxml2_source_install()
{
    source_install libxml2 "${1}" --without-python
}


function libxslt_source_install()
{
    source_install libxslt "${1}" --without-python
}


function easyrng_source_install()
{
    source_install easyrng "${1}"
}


function xmimsim_source_install()
{
    source_install xmimsim "${1}"
}


function xmimsim_post_source_install()
{
    if [[ ! -f xmimsimdata.h5 && ${ARG_SKIPLONG} == true ]]; then
        cprint "Skipping xmimsim database creation"
        return
    fi
    local prefix=${1}
    if [[ ! -f ${prefix}/bin/xmimsim-db ]]; then
        cerror "${prefix}/bin/xmimsim-db: not installed"
        return
    fi
    if [[ ! -f xmimsimdata.h5 ]]; then
        ${prefix}/bin/xmimsim-db xmimsimdata.h5
    fi
    cp -f xmimsimdata.h5 ${prefix}/share/xmimsim/xmimsimdata.h5
}


function libxml2_exists()
{
    libexists libxml2
}


function libxslt_exists()
{
    libexists libxslt
}


function easyrng_exists()
{
    libexists libeasyRNG
}


function xmimsim_exists()
{
    cmdexists xmimsim
}


function libxml2_version()
{
    libversion libxml2
}


function libxslt_version()
{
    libversion libxslt
}


function easyrng_version()
{
    libversion libeasyRNG
}


function xmimsim_version()
{
    if [[ $(xmimsim_exists) == false ]]; then
        echo 0
    else
        xmimsim --version | head -1 | grep -o -E "[0-9\.]+[0-9]"
    fi
}


function require_libxml2()
{
    require_software libxml2 ${1}
}


function require_libxslt()
{
    require_software libxslt ${1}
}


function require_easyrng()
{
    require_software easyrng ${1}
}


function require_xmimsim()
{
    require_software xmimsim ${1}
}

