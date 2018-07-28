#!/bin/bash
# 
# Install GNU libraries.
# 

SCRIPT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${SCRIPT_ROOT}/funcs.sh


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
    if [[ -z ${rversion} ]];then
        gnu_all_versions ${libname} | tail -1
        return
    fi
    
    local _version
    for i in $(gnu_all_versions ${libname}); do
        _version=${i}
        if [[ $(require_new_version ${_version} ${rversion}) == false ]]; then
            echo ${_version}
            return
        fi
    done

    echo ${_version}
}


function gnu_download()
{
    if [[ ! -f ${2}.tar.gz ]]; then
        curl -kL $(gnu_url)/${1}/${2}.tar.gz --output ${2}.tar.gz
    fi
}


function gnu_install_fromsource()
{
    local restorewd=$(pwd)
    local libname=${1}
    local rversion=${2}
    
    cprint "Download ${libname} ..."
    mkdir -p ${libname}
    cd ${libname}

    local version=$(get_local_version ${rversion})
    if [[ -z ${version} ]]; then
        require_web_essentials
        version=$(gnu_latest ${libname} ${rversion})
    fi

    local sourcedir=${libname}-${version}
    if [[ $(dryrun) == false && ! -d ${sourcedir} ]]; then
        gnu_download ${libname} ${sourcedir}
        mkdir -p ${sourcedir}
        tar -xzf ${sourcedir}.tar.gz -C ${sourcedir} --strip-components=1
        rm -f ${sourcedir}.tar.gz
    fi
    cd ${sourcedir}

    mkdir -p build
    cd build
    
    local prefix=$(project_opt)/${libname}/${version}
    local prefixstr=$(project_optstr)/${libname}/${version}
    if [[ ! -d ./bin ]]; then

        cprint "Configure ${libname} for ${prefix} ..."
        if [[ $(dryrun) == false ]]; then
            # Remove local path from LD_LIBRARY_PATH
            local keep_LD_LIBRARY_PATH=${LD_LIBRARY_PATH}
            if [[ "${LD_LIBRARY_PATH: -1}" == ":" ]];then
                LD_LIBRARY_PATH=${LD_LIBRARY_PATH::-1}
            fi
    
            mexec mkdir -p ${prefix}
            ../configure --prefix=${prefix}
            
            LD_LIBRARY_PATH=${keep_LD_LIBRARY_PATH}
        fi

        cprint "Build ${libname} ..."
        if [[ $(dryrun) == false ]]; then
            make -s -j2
        fi
    fi

    cprint "Install ${libname} in ${prefix} ..."
    if [[ $(dryrun) == false ]]; then
        if [[ ! -d ${prefix}/bin ]]; then
            mmakeinstall ${libname}-${version}
        fi

        addProfile $(project_resource) "# Installed ${libname}: ${prefixstr}"
        addBinPath ${prefix}/bin
        addBinPathProfile $(project_resource) "${prefixstr}/bin"
        addLibPath ${prefix}/lib
        addLibPathProfile $(project_resource) "${prefixstr}/lib"
    fi

    cd ${restorewd}
}


