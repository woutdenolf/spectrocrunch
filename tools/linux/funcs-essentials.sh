#!/bin/bash
#

function require_build_essentials()
{
    mapt-get make build-essential
    mapt-get make checkinstall
    mapt-get make autoconf libtool
    mapt-get make libbz2-dev zlib1g-dev
    mapt-get make openssl
}


function require_web_essentials()
{
    mapt-get wget curl
    require_web_access
}


function require_pcre()
{
    mapt-get libpcre3 libpcre3-dev
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
