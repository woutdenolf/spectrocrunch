#!/bin/bash
# 
# make helpers
# 

SCRIPT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${SCRIPT_ROOT}/funcs.sh

# ============easymake============
# Description: Execute typical configure/make/make install
#              Either directory program-version must exist
#              or file program-version.tar.gz (untar and remove)
# Usage: easymake xraylib 3.3.0 --enable-python --enable-python-integration ...
function easymake()
{
    local restorewd=$(pwd)
    local program=${1}
    shift
    local version=${1}
    shift
    
    local base=${program}-${version}
    if [[ ! -d ${base} ]]; then
        mkdir -p ${base}
        tar -xzf ${base}.tar.gz -C ${base}
        while [[ $(ls ${base} -1 | wc -l) == 1 ]];do
            local _subdir=${base}/$(ls ${base} -1)
            mv ${_subdir}/* ${base}
            rm -r ${_subdir}
        done
        rm -f ${base}.tar.gz
    fi
    cd ${base}
    mkdir -p build
    cd build

    local prefix=$(project_opt)/${program}/${version}
    local prefixstr=$(project_optstr)/${program}/${version}
    if [[ ! -e "Makefile" ]]; then
        cprint "Configure ${program} (${version}) with options:"
        cprint " --prefix=\"${prefix}\"" "$@"

        # Remove local directory from LD_LIBRARY_PATH
        local keep_LD_LIBRARY_PATH=${LD_LIBRARY_PATH}
        while [[ "${LD_LIBRARY_PATH: -1}" == ":" ]];do
            LD_LIBRARY_PATH=${LD_LIBRARY_PATH::-1}
        done
        while [[ "${LD_LIBRARY_PATH:0:1}" == ":" ]];do
            LD_LIBRARY_PATH=${LD_LIBRARY_PATH:1}
        done
        
        # Make sure destination exists
        mexec mkdir -p ${prefix}
        
        # Reconfigure is needed
        #autoreconf -i
        #if [[ -f "autogen.sh" ]]; then
        #    ../autogen.sh "$@"
        #    if [[ $? != 0 ]]; then
        #        LD_LIBRARY_PATH=${keep_LD_LIBRARY_PATH}
        #        cd ${restorewd}
        #        return
        #    fi
        #fi
        if [[ ! -f "../configure" && -f "../configure.ac" ]]; then
            cd ..
            libtoolize --force
            aclocal
            autoheader
            automake --force-missing --add-missing
            autoconf
            cd build
        fi
    
        # Configure
        #../configure --help
        #return
        ../configure --prefix="${prefix}" "$@"
        
        LD_LIBRARY_PATH=${keep_LD_LIBRARY_PATH}
        if [[ $? != 0 ]]; then
            cerror "Configuring ${program} (${version}) failed"
            cd ${restorewd}
            return
        fi
    else
        cprint "Configure ${program} (${version}): already configured."
    fi

    cprint "Build ${program} (${version}) ..."
    make -s -j2
    if [[ $? != 0 ]]; then
        cerror "Building ${program} (${version}) failed"
        cd ${restorewd}
        return
    fi
    
    cprint "Install ${program} (${version}) ..."
    mmakeinstall "${program}-${version}"
    #mmakepack ${prefix}
    #mdpkg_install *.deb ${prefix}
        
    cprint "Set environment for ${program} (${version}) ..."
    addProfile $(project_resource) "# Installed ${program}: ${prefixstr}"
    addBinPath "${prefix}/bin"
    addBinPathProfile $(project_resource) "${prefixstr}/bin"
    addLibPath "${prefix}/lib"
    addLibPathProfile $(project_resource) "${prefixstr}/lib"
    addPkgConfigPath "${prefix}/lib/pkgconfig"
    addPkgConfigPathProfile $(project_resource) "${prefixstr}/lib/pkgconfig"

    cprint "Post installation ${program} (${version}) ..."
    eval ${program}_post_source_install ${prefix}
    
    cd ${restorewd}
}


# ============source_install============
# Description: Build and install with make
# Usage: source_install "xraylib" 3.3.0  --enable-python --enable-python-integration ...
#        Needs functions: xraylib_latest, xraylib_download, xraylib_build_dependencies
function source_install()
{
    local program=${1}
    shift
    local rversion=${1}
    shift
    local restorewd=$(pwd)
    local func_latest="${program}_latest"
    local func_download="${program}_download"
    local func_deps="${program}_build_dependencies"

    mkdir -p ${program}
    cd ${program}

    cprint "Determine ${program} version to be installed ..."
    local version=$(get_local_version)
    if [[ -z ${version} || $(require_new_version ${version} ${rversion}) ]]; then
        require_web_essentials
        version=$($func_latest ${rversion})
    fi
    
    local base=${program}-${version}
    echo ${base}
    if [[ $(dryrun) == false && ! -d ${base} ]]; then
        cprint "Download ${program} ..."
        eval $func_download ${base}
    fi

    if [[ $(dryrun) == false ]]; then
        cd ${restorewd}
        eval $func_deps
        cd ${program}
        easymake ${program} ${version} "$@"
    fi

    cd ${restorewd}
}


# ============require_software============
# Description: Get the latest version
# Usage: require_software "swig" 3
#        Needs functions: swig_version
#                         swig_exists
#                         swig_system_install
#                         swig_source_install
#                         swig_run_dependencies
#                         swig_build_dependencies
#                         swig_latest
#                         swig_download
function require_software()
{
    local program=${1}
    local rversion=${2}
    
    cprintstart
    cprint "Verify ${program} ${rversion} ..."

    # Check version
    if [[ $(require_new_version $(${program}_version) ${rversion}) == false ]]; then
        cprint "${program} version $(${program}_version) will be used"
        cprintend
        return
    fi
    
    # Try system installation
    if [[ $(dryrun) == false ]]; then
        cprint "System install ${program} ${rversion} ..."
        eval ${program}_system_install ${rversion}

        # Check version
        if [[ $(require_new_version $(${program}_version) ${rversion}) == false ]]; then
            cprint "${program} version $(${program}_version) will be used"
            cprintend
            return
        fi
    fi
    
    # Run dependencies
    if [[ $(dryrun) == false ]]; then
        cprint "Install run dependencies of ${program} ${rversion} ..."
        eval ${program}_run_dependencies
    fi
    
    # Install from source
    cprint "Source install ${program} ${rversion} ..."
    eval ${program}_source_install ${rversion}

    # Check version
    if [[ $(require_new_version $(${program}_version) ${rversion}) == false ]]; then
        cprint "${program} version $(${program}_version) will be used"
    else
        if [[ $(${program}_exists) == false ]]; then
            cerror "${program} is not installed"
        else
            cerror "${program} version $(${program}_version) will be used but ${rversion} is required"
        fi
    fi

    cprintend
}


# ============latest_version============
# Description: Get the latest version
# Usage: latest_version "swig_all_versions" 3
function latest_version()
{
    local _all_versions=${1}
    local rversion=${2}

    # Last version
    if [[ -z ${rversion} ]];then
        ($_all_versions) | tail -1
        return
    fi
    
    # Last version with the same base
    local _version
    for i in $($_all_versions); do
        if [[ $(require_new_version ${i} ${rversion}) == false ]]; then
            if [[ ${i} == ${rversion} || ${i} =~ ${rversion}[^0-9] ]]; then
                _version=${i}
            fi
        fi
    done

    # Last version equal or better
    if [[ -z ${_version} ]];then
        for i in $($_all_versions); do
            if [[ $(require_new_version ${i} ${rversion}) == false ]]; then
                _version=${i}
            fi
        done
    fi

    echo ${_version}
}


# ============versions_from_site============
# Description: Get all releases from site
# Usage: versions_from_site "ftp://xmlsoft.org/libxslt/" "libxslt-[0-9\.]+[0-9]\.tar\.gz$"
function versions_from_site()
{
    curl -sL ${1} | grep -E -o ${2} | grep -E -o "[0-9\.]+[0-9]" | sort --version-sort
}


# ============versions_from_github============
# Description: Get all releases on github
# Usage: versions_from_github "tschoonj" "easyRNG" "easyRNG-[0-9\.]+"
function versions_from_github()
{
    local header=""
    if [[ ! -z ${GITHUB_TOKEN} ]];then
        header="Authorization: token ${GITHUB_TOKEN}"
    fi
    curl -H "${header}" --silent "https://api.github.com/repos/${1}/${2}/tags" | grep -o -E "\"name\": \"${3}\"" | grep -E -o "[0-9\.]+[0-9]" | sort --version-sort
}

