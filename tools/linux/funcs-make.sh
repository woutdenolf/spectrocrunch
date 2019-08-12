#!/bin/bash
# 
# make helpers
# 

SCRIPT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${SCRIPT_ROOT}/funcs.sh


# ============require_new_version============
# Description: a new version is required when current < required
# Usage: [[ $(require_new_version currentversion requiredversion) ]]
function require_new_version()
{
    local currentv=${1}
    local requiredv=${2}

    if [[ ${currentv} == 0 || -z ${currentv} ]]; then
        echo true # does not exist
        return
    fi

    if [[ -z ${requiredv} ]]; then
        echo false # no specific version required
        return
    fi

    dpkg --compare-versions ${currentv} "lt" ${requiredv}
    if [[ $? == 0 ]]; then
        echo true # current < required
    else
        echo false # current >= required
    fi
}


# ============require_new_version_strict============
# Description: a new version is required when current != required (common depth)
# Usage: [[ $(require_new_version_strict currentversion requiredversion) ]]
function require_new_version_strict()
{
    local currentv=${1}
    local requiredv=${2}

    if [[ ${currentv} == 0 ]]; then
        echo true # does not exist
        return
    fi

    if [[ -z ${requiredv} ]]; then
        echo false # no specific version required
        return
    fi

    # Same version depth
    local n=$(nchar_occurence ${currentv} .)
    local nreq=$(nchar_occurence ${requiredv} .)
    if [[ ${nreq} -lt ${n} ]];then
        n=${nreq}
    fi
    n=$((n+1))
    currentv=$(echo ${currentv} | grep -o -E "[0-9]+" | head -${n})
    requiredv=$(echo ${requiredv} | grep -o -E "[0-9]+" | head -${n})
    
    currentv=$(echo ${currentv})
    requiredv=$(echo ${requiredv})
    currentv=${currentv// /.}
    requiredv=${requiredv// /.}

    # Compare
    dpkg --compare-versions "${currentv}" "eq" "${requiredv}"
    if [[ $? == 0 ]]; then
        echo false # current == required
    else
        echo true # current != required
    fi
}


# ============get_local_version============
# Description: 
# Usage: version=$(get_local_version requiredversion)
function get_local_version()
{
    local requiredv=${1}

    for dirname in $(ls -d */ | sort -V); do
        local version=$(echo ${dirname} | grep -E -o "[0-9\.]+[0-9]")
        if [[ $(require_new_version ${version} ${requiredv}) == false ]]; then
            echo ${version}
            return
        fi
    done
}


# ============get_local_version_strict============
# Description: 
# Usage: version=$(get_local_version_strict requiredversion)
function get_local_version_strict()
{
    local requiredv=${1}

    for dirname in $(ls -d */ | sort -V); do
        local version=$(echo ${dirname} | grep -E -o "[0-9\.]+[0-9]")
        if [[ $(require_new_version_strict ${version} ${requiredv}) == false ]]; then
            echo ${version}
            return
        fi
    done
}


# ============easymake_prefix============
# Description: 
# Usage: easymake_prefix xraylib 3.3.0
function easymake_prefix()
{
    local program=${1}
    local version=${2}
    echo $(project_opt)/${program}/${version}
}


# ============easymake_prefixstr============
# Description: 
# Usage: easymake_prefixstr xraylib 3.3.0
function easymake_prefixstr()
{
    local program=${1}
    local version=${2}
    echo $(project_optstr)/${program}/${version}
}


# ============mmakeinstall============
# Description: Execute make install
# Usage: mmakeinstall pkgname-version
function mmakeinstall()
{
    local name=${1}
    if [[ -z ${name} ]];then
        name=$(randomstring 6)
    fi
    if [[ $(install_systemwide) == true ]]; then
        if [[ $(cmdexists "checkinstall") == true ]]; then    
            sudo -E checkinstall -y --pkgname "${name}-checkinstall"
            dpkg -l ${name}-checkinstall
            echo "Remove with \"dpkg -r ${name}-checkinstall\""
        else
            sudo -E make install -s
        fi
    else
        make install -s
    fi
}


# ============mmakepack============
# Description: Execute make install
# Usage: mmakepack pkgname-version
function mmakepack()
{
    local name=${1}
    if [[ -z ${name} ]];then
        name=$(randomstring 6)
    fi
    checkinstall -y --pkgname "${name}-checkinstall" --install=no
}


# ============mdpkg_install============
# Description: dpkg without prompt
# Usage: mdpkg_install package.deb ${prefix}
function mdpkg_install()
{
    local package="$1"
    local prefix="$2"
    local extension=${package: -4}
    if [[ ${extension} == ".deb" ]]; then
        if [[ $(install_systemwide) == true ]]; then
            sudo -E dpkg -i ${package}
        else
            dpkg -x ${package} ${prefix}
        fi
    elif [[ ${extension} == ".rpm" ]]; then
        if [[ $(install_systemwide) == true ]]; then
            sudo -E rpm -i ${package}
        else
            local restorewd=$(pwd)
            cd ${prefix}
            rpm2cpio "${restorewd}/${package}" | cpio -id
            cd ${restorewd}
        fi
    else
        echo "Skip ${package} installation (unknown package extension)"
    fi
}


# ============cprint_makeenv============
# Description: 
# Usage: 
function cprint_makeenv()
{
    cprint "Binaries: PATH=${PATH}"
    cprint "Headers: CPATH=${CPATH}"
    cprint "Shared libraries: LD_LIBRARY_PATH=${LD_LIBRARY_PATH}"
    cprint "Static libraries: LIBRARY_PATH=${LIBRARY_PATH}"
    cprint "Library metadata: PKG_CONFIG_PATH=${PKG_CONFIG_PATH}"
}


# ============removelocal_path============
# Description: 
# Usage: 
function removelocal_path()
{
    local _path=${1}
    while [[ "${_path: -1}" == ":" ]];do
        _path=${_path::-1}
    done
    while [[ "${_path:0:1}" == ":" ]];do
        _path=${_path:1}
    done
    echo ${_path}
}


# ============makeenv_removelocal============
# Description: 
# Usage: 
function makeenv_removelocal()
{
    export PATH=$(removelocal_path ${PATH})
    export CPATH=$(removelocal_path ${CPATH})
    export LD_LIBRARY_PATH=$(removelocal_path ${LD_LIBRARY_PATH})
    export LIBRARY_PATH=$(removelocal_path ${LIBRARY_PATH})
    export PKG_CONFIG_PATH=$(removelocal_path ${PKG_CONFIG_PATH})
}


# ============easymake_makeenv============
# Description: 
# Usage: 
function easymake_makeenv()
{
    local program=${1}
    local version=${2}
    local prefix=$(easymake_prefix ${program} ${version})
    local prefixstr=$(easymake_prefixstr ${program} ${version})

    addProfile $(project_resource) "# Installed ${program}: ${prefixstr}"
    if [[ -d "${prefix}/bin/${program}" ]];then
        addBinPath "${prefix}/bin/${program}"
        addBinPathProfile $(project_resource) "${prefixstr}/bin/${program}"
    else
        addBinPath "${prefix}/bin"
        addBinPathProfile $(project_resource) "${prefixstr}/bin"
    fi
    if [[ -d "${prefix}/lib/${program}" ]];then
        addLibPath "${prefix}/lib/${program}"
        addLibPathProfile $(project_resource) "${prefixstr}/lib/${program}"
    else
        addLibPath "${prefix}/lib"
        addLibPathProfile $(project_resource) "${prefixstr}/lib"
    fi
    if [[ -d "${prefix}/include/${program}" ]];then
        addInclPath "${prefix}/include/${program}"
        addInclPathProfile $(project_resource) "${prefixstr}/include/${program}"
    else
        addInclPath "${prefix}/include"
        addInclPathProfile $(project_resource) "${prefixstr}/include"
    fi
    if [[ -d "${prefix}/lib/${program}/pkgconfig" ]];then
        addPkgConfigPath "${prefix}/lib/${program}/pkgconfig"
        addPkgConfigPathProfile $(project_resource) "${prefixstr}/lib/${program}/pkgconfig"
    else
        addPkgConfigPath "${prefix}/lib/pkgconfig"
        addPkgConfigPathProfile $(project_resource) "${prefixstr}/lib/pkgconfig"
    fi
}


# ============easymake_configure============
# Description: 
# Usage: 
function easymake_configure()
{
    local program=${1}
    local version=${2}
    local cfgparams="${@:3}"
    local prefix=$(easymake_prefix ${program} ${version})

    if [[ -e "Makefile" ]]; then
        cprint "Configure ${program} (${version}): already configured."
        cprint " --prefix=\"${prefix}\" ${cfgparams}"
    else
        cprint "Configure ${program} (${version}) with options:"
        cprint " --prefix=\"${prefix}\" ${cfgparams}" 
        
        # In case configure file does not exist
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
        ../configure --help
        ../configure --prefix="${prefix}" ${cfgparams}
        
        if [[ $? != 0 ]]; then
            cerror "Configuring ${program} (${version}) failed"
            LD_LIBRARY_PATH=${keep_LD_LIBRARY_PATH}
            cd ${restorewd}
            return
        fi
    fi
}


# ============easymake_build============
# Description: 
# Usage: 
function easymake_build()
{
    local program=${1}
    local version=${2}
    make -s -j$(($(nproc)+1))
    return $?
}


# ============easymake_install============
# Description: 
# Usage: 
function easymake_install()
{
    local program=${1}
    local version=${2}
    local prefix=$(easymake_prefix ${program} ${version})
    mmakeinstall "${program}-${version}"
    #mmakepack ${prefix}
    #mdpkg_install *.deb ${prefix}
}


# ============easymake============
# Description: Execute typical configure/make/make install
#              Either directory program-version must exist
#              or file program-version.tar.gz (untar and remove)
# Usage: easymake xraylib 3.3.0 --enable-python --enable-python-integration ...
function easymake()
{
    local restorewd=$(pwd)
    local program=${1}
    local version=${2}
    local cfgparams="${@:3}"
    local func_configure="${program}_configure"
    local func_build="${program}_build"
    local func_install="${program}_install"
    local func_post="${program}_post"
    local prefix=$(easymake_prefix ${program} ${version})
    cfgparams=$(eval echo ${cfgparams})

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

    # Make sure destination exists
    mexec mkdir -p ${prefix}
    
    # Remove local directory env paths
    makeenv_removelocal
    cprint_makeenv

    cprint "Configure ${program} (${version}) with options:"
    local prefix=$(easymake_prefix ${program} ${version})
    cfgparams=$(eval echo ${cfgparams})
    cprint "${cfgparams}"
    if [[ $(cmdexists ${func_configure}) == true ]];then
        eval ${func_configure} ${program} ${version} ${cfgparams}
    else
        easymake_configure ${program} ${version} ${cfgparams}
    fi

    cprint "Build ${program} (${version}) ..."
    if [[ $(cmdexists ${func_build}) == true ]];then
        eval ${func_build} ${program} ${version}
    else
        easymake_build ${program} ${version}
    fi
    if [[ $? != 0 ]]; then
        cerror "Building ${program} (${version}) failed"
        cd ${restorewd}
        return
    fi
    
    cprint "Install ${program} (${version}) ..."
    if [[ $(cmdexists ${func_install}) == true ]];then
        eval ${func_install} ${program} ${version}
    else
        easymake_install ${program} ${version}
    fi

    cprint "Set environment for ${program} (${version}) ..."
    if [[ $(cmdexists ${func_makeenv}) == true ]];then
        eval ${func_makeenv} ${program} ${version}
    else
        easymake_makeenv ${program} ${version}
    fi

    cprint "Post installation ${program} (${version}) ..."
    if [[ $(cmdexists ${func_post}) == true ]];then
        eval ${func_post} ${program} ${version}
    fi
    
    cd ${restorewd}
}


# ============source_install============
# Description: Build and install with make
# Usage: source_install "xraylib" 3.3.0  --enable-python --enable-python-integration ...
#        Needs functions: xraylib_configure (optional)
#                         xraylib_latest
#                         xraylib_download
#                         xraylib_build_dependencies (optional)
function source_install()
{
    local program=${1}
    local rversion=${2}
    local restorewd=$(pwd)
    local func_latest="${program}_latest"
    local func_download="${program}_download"
    local func_deps="${program}_build_dependencies"
    local cfgparams="${@:3}"

    mkdir -p ${program}
    cd ${program}

    cprint "Determine ${program} version to be installed ..."
    local version=$(get_local_version)
    if [[ -z ${version} || $(require_new_version ${version} ${rversion}) ]]; then
        require_web_essentials
        version=$($func_latest ${rversion})
    fi
    
    if [[ -z ${version} ]];then
        version="master"
    fi

    local base=${program}-${version}
    if [[ $(dryrun) == false && ! -d ${base} ]]; then
        cprint "Download ${program} (${version}) ..."
        eval $func_download ${base}
    fi

    if [[ $(dryrun) == false ]]; then
        cd ${restorewd}
        eval $func_deps
        cd ${program}
        easymake ${program} ${version} ${cfgparams}
    fi

    cd ${restorewd}
}


# ============require_software============
# Description: Get the latest version
# Usage: require_software "swig" 3
#        Needs functions: swig_version
#                         swig_exists
#                         swig_system_install (optional)
#                         swig_source_install
#                         swig_run_dependencies (optional)
#                         swig_build_dependencies (optional)
#                         swig_configure (optional)
#                         swig_build (optional)
#                         swig_install (optional)
#                         swig_environment (optional)
#                         swig_post (optional)
#                         swig_latest
#                         swig_download
function require_software()
{
    local program=${1}
    local rversion=${2}
    cprintstart "Require ${program} ${rversion}"

    # Check version
    if [[ $(require_new_version $(${program}_version) ${rversion}) == false ]]; then
        cprint "${program} version $(${program}_version) will be used"
        cprintend "Require ${program} ${rversion}"
        return
    fi
    
    if [[ $(dryrun) == false ]]; then
        # Try system installation
        cprint "System install ${program} ${rversion} ..."
        eval ${program}_system_install ${rversion}

        # Check version
        if [[ $(require_new_version $(${program}_version) ${rversion}) == false ]]; then
            cprint "${program} version $(${program}_version) will be used"
            cprintend "Require ${program} ${rversion}"
            return
        fi

        # Run dependencies
        cprint "Install run dependencies of ${program} ${rversion} ..."
        eval ${program}_run_dependencies
    
        # Install from source
        cprint "Source install ${program} ${rversion} ..."
        eval ${program}_source_install ${rversion}
    fi

    # Check version
    if [[ $(require_new_version $(${program}_version) ${rversion}) == false ]]; then
        cprint "${program} version $(${program}_version) will be used"
    else
        local _exists
        if [[ $(cmdexists "${program}_exists") == true ]];then
            _exists=$(${program}_exists)
        else
            _exists=$(cmdexists ${program})
        fi
        if [[ ${_exists} == false ]]; then
            cerror "${program} is not installed"
        else
            cerror "${program} version $(${program}_version) will be used but ${rversion} is required"
        fi
    fi

    cprintend "Require ${program} ${rversion}"
}


# ============latest_version============
# Description: Get the latest version
# Usage: latest_version "swig_all_versions" 3
function latest_version()
{
    local lst=($(eval ${1}))
    local rversion=${2}

    # Last version when no version requested
    if [[ -z ${rversion} ]];then
        echo ${lst[-1]}
        return
    fi
    
    # Last version with the same base
    local _version=""
    for i in ${lst[@]}; do
        if [[ $(require_new_version ${i} ${rversion}) == false ]]; then
            if [[ ${i} == ${rversion} || ${i} =~ ${rversion}[^0-9] ]]; then
                _version=${i}
            fi
        fi
    done

    # Last version equal or better
    if [[ -z ${_version} ]];then
        for i in ${lst[@]}; do
            if [[ $(require_new_version ${i} ${rversion}) == false ]]; then
                _version=${i}
            fi
        done
    fi

    # Last version 
    if [[ -z ${_version} ]];then
        _version=${lst[-1]}
    fi

    echo ${_version}
}


# ============versions_from_site============
# Description: Get all releases from site
# Usage: versions_from_site "ftp://xmlsoft.org/libxslt" "libxslt-[0-9\.]+[0-9]\.tar\.gz$"
function versions_from_site()
{
    curl -slL ${1} | grep -E -o ${2} | grep -E -o "[0-9\.]+[0-9]" | sort --version-sort
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

