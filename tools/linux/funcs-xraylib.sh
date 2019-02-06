#!/bin/bash
# 
# Install xraylib.
# 

SCRIPT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${SCRIPT_ROOT}/funcs.sh
source ${SCRIPT_ROOT}/funcs-swig.sh
source ${SCRIPT_ROOT}/funcs-python.sh


function xraylib_build_dependencies()
{
    local tmp=$(pwd)
    cd ${1}

    require_build_essentials
    require_pythondev
    require_swig 3
    pip_install cython

    cd ${tmp}
}


function xraylib_latest()
{
    local link=""
    if [[ -z ${GITHUB_TOKEN} ]];then
        link="https://api.github.com/repos/tschoonj/xraylib/tags?access_token=${GITHUB_TOKEN}"
    else
        link="https://api.github.com/repos/tschoonj/xraylib/tags"
    fi
    curl --silent "https://api.github.com/repos/tschoonj/xraylib/tags" | grep -o -E "xraylib-[0-9\.]+" | head -1 | grep -E -o "[0-9\.]+"
}


function xraylib_download()
{
    if [[ ! -f ${1}.tar.gz ]]; then
        curl -L https://github.com/tschoonj/xraylib/archive/${1}.tar.gz --output ${1}.tar.gz
    fi
}


function xraylib_install_fromsource()
{
    #if [[ ! -d xraylib && ${ARG_SKIPLONG} == true ]]; then
    #    cprint "Skipping xraylib installation"
    #    return
    #fi

    local restorewd=$(pwd)

    cprint "Download xraylib ..."
    mkdir -p xraylib
    cd xraylib

    local version=$(get_local_version)
    if [[ -z ${version} ]]; then
        require_web_essentials
        version=$(xraylib_latest)
    fi
    
    local sourcedir=xraylib-${version}
    if [[ $(dryrun) == false && ! -d ${sourcedir} ]]; then
        xraylib_download ${sourcedir}
        mkdir -p ${sourcedir}
        tar -xzf ${sourcedir}.tar.gz -C ${sourcedir} --strip-components=1
        rm -f ${sourcedir}.tar.gz
    fi
    cd ${sourcedir}

    local prefix=$(project_opt)/xraylib/${version}
    local prefixstr=$(project_optstr)/xraylib/${version}
    if [[ ! -f python/.libs/_xraylib.so ]]; then

        cprint "Configure xraylib for ${prefix} ..."
        if [[ $(dryrun) == false ]]; then
            xraylib_build_dependencies ${restorewd}

            mexec mkdir -p ${prefix}
            autoreconf -i
            ./configure --prefix="${prefix}" \
                        --enable-python \
                        --enable-python-integration \
                        --disable-java \
                        --disable-lua \
                        --disable-ruby \
                        --disable-php \
                        --disable-pascal \
                        --disable-idl \
                        --disable-perl \
                        --disable-fortran2003 \
                        PYTHON="$(python_full_bin)"
            $(pip_bin) freeze | grep numpy > requirements.txt
        fi

        cprint "Build xraylib ..."
        if [[ $(dryrun) == false ]]; then
            make -s -j2
        fi
    fi

    cprint "Install xraylib in ${prefix} ..."
    if [[ $(dryrun) == false ]]; then
        pip_install -r requirements.txt
        mmakeinstall xraylib-${version}
        #mmakepack xraylib-${version}
        #mdpkg_install *.deb ${prefix}

        addProfile $(project_resource) "# Installed xraylib: ${prefixstr}"
        addBinPath ${prefix}/bin
        addBinPathProfile $(project_resource) "${prefixstr}/bin"
        addLibPath ${prefix}/lib
        addLibPathProfile $(project_resource) "${prefixstr}/lib"
    fi

    cd ${restorewd}
}


function require_xraylib()
{
    cprintstart
    cprint "Verify xraylib ..."

    # Requirements (for running)
    require_python

    # Check
    if [[ $(python_hasmodule xraylib) == true ]]; then
        cprint "Python module \"xraylib\" is installed"
        cprintend
        return
    fi

    # Install from source
    xraylib_install_fromsource

    # Check
    if [[ $(python_hasmodule xraylib) == true ]]; then
        cprint "Python module \"xraylib\" is installed"
    else
        cprint "Python module \"xraylib\" is NOT installed"
    fi

    cprintend
}


