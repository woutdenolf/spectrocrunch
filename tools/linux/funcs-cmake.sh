#!/bin/bash
# 
# Install cmake.
# 

SCRIPT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${SCRIPT_ROOT}/funcs.sh


function cmake_url()
{
    echo "https://cmake.org/files/LatestRelease/"
}


function cmake_latest()
{
    curl -sL $(cmake_url) | grep -E -o "cmake-[0-9\.]+[0-9]" | tail -1 | grep -E -o "[0-9\.]+[0-9]"
}


function cmake_download()
{
    if [[ ! -f ${1}.tar.gz ]]; then
        curl -L $(cmake_url)/${1}.tar.gz --output ${1}.tar.gz
    fi
}


function cmake_build_dependencies()
{
    require_build_essentials
}


function cmake_install_fromsource()
{
    local restorewd=$(pwd)

    cprint "Download cmake ..."
    mkdir -p cmake
    cd cmake

    local version=$(get_local_version ${1})
    if [[ -z ${version} ]]; then
        require_web_essentials
        local version=$(cmake_latest)
    fi

    local sourcedir=cmake-${version}
    if [[ $(dryrun) == false && ! -d ${sourcedir} ]]; then
        cmake_download ${sourcedir}
        mkdir -p ${sourcedir}
        tar -xzf ${sourcedir}.tar.gz -C ${sourcedir} --strip-components=1
        rm -f ${sourcedir}.tar.gz
    fi
    cd ${sourcedir}

    local prefix=$(project_opt)/cmake/${version}
    if [[ ! -f ./bin/cmake ]]; then

        cprint "Configure cmake for ${prefix} ..."
        if [[ $(dryrun) == false ]]; then
            cmake_build_dependencies

            mexec mkdir -p ${prefix}
            ./configure --prefix=${prefix}
        fi

        cprint "Build cmake ..."
        if [[ $(dryrun) == false ]]; then
            make -s -j2
        fi
    fi

    cprint "Install cmake in ${prefix} ..."
    if [[ $(dryrun) == false ]]; then
        if [[ ! -f ${prefix}/bin/cmake ]]; then
            mmakeinstall cmake-${version}
        fi

        # Add path just for this session
        addBinPath ${prefix}/bin
        #addProfile $(project_resource) "# Installed cmake: ${prefixstr}"
        #addBinPath ${prefix}/bin
        #addBinPathProfile $(project_resource) "${prefixstr}/bin"

    fi

    cd ${restorewd}
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
    cprintstart
    cprint "Verify cmake ${1} ..."

    # Try system installation
    if [[ $(cmdexists cmake) == false ]]; then
        mapt-get install cmake
    fi

    # Check version
    if [[ $(require_new_version $(cmake_version) ${1}) == false ]]; then
        cprint "cmake version $(cmake_version) will be used"
        cprintend
        return
    fi

    # Install from source
    cmake_install_fromsource ${1}

    # Check version
    if [[ $(require_new_version $(cmake_version) ${1}) == false ]]; then
        cprint "cmake version $(cmake_version) will be used"
    else
        if [[ $(cmdexists cmake) == false ]]; then
            cerror "cmake is not installed"
        else
            cerror "cmake version $(cmake_version) will be used but ${1} is required"
        fi
    fi

    cprintend
}



