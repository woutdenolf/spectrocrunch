#!/bin/bash
# 
# Install swig.
# 

SCRIPT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${SCRIPT_ROOT}/funcs.sh


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
    if [[ -z ${1} ]];then
        swig_all_version | tail -1
        return
    fi
    
    local _version
    for i in $(swig_all_versions); do
        _version=${i}
        if [[ $(require_new_version ${_version} ${1}) == false ]]; then
            echo ${_version}
            return
        fi
    done

    echo ${_version}
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


function swig_install_fromsource()
{
    local restorewd=$(pwd)

    cprint "Download swig ..."
    mkdir -p swig
    cd swig

    local version=$(get_local_version ${1})
    if [[ -z ${version} ]]; then
        require_web_essentials
        version=$(swig_latest ${1})
    fi

    local sourcedir=swig-${version}
    if [[ $(dryrun) == false && ! -d ${sourcedir} ]]; then
        swig_download ${sourcedir}
        mkdir -p ${sourcedir}
        tar -xzf ${sourcedir}.tar.gz -C ${sourcedir} --strip-components=1
        rm -f ${sourcedir}.tar.gz
    fi
    cd ${sourcedir}

    local prefix=$(project_opt)/swig/${version}
    if [[ ! -f ./bin/swig ]]; then

        cprint "Configure swig for ${prefix} ..."
        if [[ $(dryrun) == false ]]; then
            swig_build_dependencies

            mexec mkdir -p ${prefix}
            ./configure --prefix=${prefix}
        fi

        cprint "Build swig ..."
        if [[ $(dryrun) == false ]]; then
            make -s -j2
        fi
    fi

    cprint "Install swig in ${prefix} ..."
    if [[ $(dryrun) == false ]]; then
        if [[ ! -f ${prefix}/bin/swig ]]; then
            mmakeinstall swig-${version}
        fi

        # Add path just for this session
        addBinPath ${prefix}/bin
        #addProfile $(project_resource) "# Installed swig: ${prefixstr}"
        #addBinPath ${prefix}/bin
        #addBinPathProfile $(project_resource) "${prefixstr}/bin"
    fi

    cd ${restorewd}
}


function swig_version()
{
    if [[ $(cmdexists swig) == false ]]; then
        echo 0
    else
        swig -version | head -2 | tail -1 | awk '{print $3}'
    fi
}


function require_swig()
{
    cprintstart
    cprint "Verify swig ${1} ..."
    hash -d swig # swig2/swig3 conflicts

    # Try system installation
    if [[ $(cmdexists swig) == false ]]; then
        mapt-get install swig
    fi

    # Check version
    if [[ $(require_new_version $(swig_version) ${1}) == false ]]; then
        cprint "swig version $(swig_version) will be used"
        cprintend
        return
    fi

    # Install from source
    swig_install_fromsource ${1}

    # Check version
    if [[ $(require_new_version $(swig_version) ${1}) == false ]]; then
        cprint "swig version $(swig_version) will be used"
    else
        if [[ $(cmdexists swig) == false ]]; then
            cerror "swig is not installed"
        else
            cerror "swig version $(swig_version) will be used but ${1} is required"
        fi
    fi

    cprintend
}


