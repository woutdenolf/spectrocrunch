#!/bin/bash
#
# Essentials for installation and downloading
#

function require_build_essentials()
{
    mapt-get install build-essential
    mapt-get install checkinstall
    mapt-get install autoconf
    mapt-get install libtool
    mapt-get install pkgconf
}


function require_web_essentials()
{
    mapt-get install wget curl
    require_web_access
}


function require_openssl()
{
    mapt-get install openssl
    mapt-get install openssl-dev
    mapt-get install openssl-devel
    mapt-get install libssl-dev
}


function require_web_access()
{
    if [[ "$(dnsdomainname)" == "esrf.fr" ]]; then
        if [[ -z $http_proxy ]]; then
            export http_proxy="http://proxy.esrf.fr:3128"
        fi
        if [[ -z $http_proxy ]]; then
            export https_proxy="http://proxy.esrf.fr:3128"
        fi
        export use_proxy=yes
    fi
}
