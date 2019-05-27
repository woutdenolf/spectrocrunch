#!/bin/bash
# 
# Install cmake.
# 

SCRIPT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${SCRIPT_ROOT}/funcs-make.sh


function cmake_url()
{
    echo "https://cmake.org/files"
}


function cmake_all_versions()
{
    curl -sL $(cmake_url) | grep -E -o ">v[0-9\.]+[0-9]" | grep -E -o "[0-9\.]+[0-9]" | sort --version-sort
}


function cmake_all_subversions()
{
    curl -sL $(cmake_url)/v${1} | grep -E -o ">cmake-[0-9\.]+.tar.gz<" | grep -E -o "[0-9\.]+[0-9]" | sort --version-sort
}


function cmake_extractmainv()
{
    echo ${1} | grep -E -o "[0-9]+(\.[0-9]+)?" | head -1
}


function cmake_latest()
{
    local rversion=${1}
    local lst=($(cmake_all_versions))

    # Last version when no version requested
    if [[ -z ${rversion} ]];then
        echo ${lst[-1]}
        return
    fi
    local mversion=$(cmake_extractmainv ${rversion})
    local lst2
    for i in ${lst[@]}; do
        if [[ $(require_new_version ${i} ${mversion}) == false ]]; then
            lst2=($(cmake_all_subversions ${i}))
            if [[ ${mversion} == ${rversion} ]];then
                lst2=(${lst2[-1]})  
            fi
            for j in ${lst2[@]}; do
                if [[ $(require_new_version ${j} ${rversion}) == false ]]; then
                    echo ${j}
                    return
                fi
            done
        fi
    done

    cmake_all_subversions ${lst[-1]}
}


function cmake_download()
{
    if [[ ! -f ${1}.tar.gz ]]; then
        local _mversion=$(cmake_extractmainv ${1})
        curl -L $(cmake_url)/v${_mversion}/${1}.tar.gz --output ${1}.tar.gz
    fi
}


function cmake_build_dependencies()
{
    require_build_essentials
    require_openssl
    mapt-get install libncurses5-dev
    mapt-get install qt4-qmake
}


function cmake_system_install()
{
    mapt-get install cmake
}


function cmake_configure()
{
    ../bootstrap "$@"
}


function cmake_source_install()
{
    source_install cmake "${1}" -- -DCMAKE_USE_OPENSSL:BOOL=ON
}


function cmake_version()
{
    if [[ $(cmdexists cmake) == false ]]; then
        echo 0
    else
        cmake --version | head -1 | awk '{print $3}'
    fi
}


function require_cmake()
{
    require_software cmake ${1}
}



