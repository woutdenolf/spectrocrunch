#!/bin/bash
# 
# Install swig.
# 

SCRIPT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $SCRIPT_ROOT/funcs.sh


function swig_url()
{
    echo "https://sourceforge.net/projects/swig/files/swig/"
}


function swig_latest()
{
    curl -sL $(swig_url) | grep -E -o "swig-[0-9\.]+" | head -1 | grep -E -o "[0-9\.]+[0-9]"
}


function swig_download()
{
    if [[ ! -f ${1}.tar.gz ]]; then
        curl -L $(swig_url)/${1}/${1}.tar.gz --output ${1}.tar.gz
    fi
}


function swig_build_dependencies()
{
    require_build_essentials
    mapt-get "install libpcre3-dev"
}


function swig_install_fromsource()
{
    local restorewd=$(pwd)

    cprint "Download swig ..."
    mkdir -p swig
    cd swig

    require_web_essentials
    local version=$(swig_latest)
    local sourcedir=swig-${version}
    if [[ $(dryrun) == false && ! -d ${sourcedir} ]]; then
        swig_download ${sourcedir}
        mkdir -p ${sourcedir}
        tar -xzf ${sourcedir}.tar.gz -C ${sourcedir} --strip-components=1
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
        mmakeinstall swig-${version}

        # Add path just for this session
        addBinPath ${prefix}/bin
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
    hash -d swig

    # Try system installation
    if [[ $(cmdexists swig) == false ]]; then
        mapt-get "install swig"
    fi

    # Check version
    if [[ $(require_new_version $(swig_version) ${1}) == false ]]; then
        cprint "swig version $(swig_version) will be used"
        cprintend
        return
    fi

    # Install from source
    swig_install_fromsource

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


