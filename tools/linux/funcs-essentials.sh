#!/bin/bash
#

function require_build_essentials()
{
    mapt-get install make build-essential
    mapt-get install make checkinstall
    mapt-get install make autoconf libtool
    mapt-get install make libbz2-dev zlib1g-dev
    mapt-get install make openssl openssl-devel
}


function require_web_essentials()
{
    mapt-get install wget curl
    require_web_access
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
