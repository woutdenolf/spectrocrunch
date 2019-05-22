#!/bin/bash
# 
# Helper functions.
# 

SCRIPT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${SCRIPT_ROOT}/funcs-string.sh
source ${SCRIPT_ROOT}/funcs-essentials.sh

# ============fullpath============
# Description: expand path
function fullpath()
{
    eval echo ${1}
}


# ============cprint============
# Description: output to stdout in color
# Usage: cprint "..."
function cprint()
{
    local hcol='\033[0;36m'
    local ncol='\033[0m'
    echo -e "${hcol}$@ ${ncol}"
}


# ============cerror============
# Description: output to stdout in color
# Usage: cerror "..."
function cerror()
{
    local hcol='\033[0;31m'
    local ncol='\033[0m'
    echo -e "${hcol}$@ ${ncol}"
}


# ============cprintstart============
# Description: output to stdout in color
# Usage: cprintstart
function cprintstart()
{
    echo ""
    echo ""
    echo ""
    cprint "======================"
}


# ============cprintend============
# Description: output to stdout in color
# Usage: cprintend
function cprintend()
{
    cprint "======================"
    echo ""
    echo ""
    echo ""
}


# ============system_privileges============
# Description: check for root access
# Usage: [[ $(system_privileges) ]]
function system_privileges()
{
    if [[ -z "$((sudo -n true) 2>&1)" ]]; then
        echo true 
    else
        echo false
    fi
}


# ============install_systemwide============
# Description: check for user or system installation
# Usage: [[ $(install_systemwide) == true ]]
#        install_systemwide reset true
function install_systemwide()
{
    if [[ -z ${INSTALL_SYSTEMWIDE} || "$1" == "reset" ]]; then
        if [[ -z ${2} ]]; then
            INSTALL_SYSTEMWIDE=$(system_privileges)
        else
            INSTALL_SYSTEMWIDE=${2}
        fi
        if [[ "$1" == "reset" ]];then
            return
        fi
    fi

    if [[ $(system_privileges) == false ]]; then
        echo false
    else
        echo ${INSTALL_SYSTEMWIDE}
    fi
}


# ============dryrun============
# Description: check whether we are doing a dry run
# Usage: [[ $(dryrun) == true ]]
#        dryrun reset true
function dryrun()
{
    if [[ -z ${DRYRUN} || "$1" == "reset" ]]; then
        if [[ -z ${2} ]]; then
            DRYRUN=true
        else
            DRYRUN=${2}
        fi
        if [[ "$1" == "reset" ]];then
            return
        fi
    fi

    echo ${DRYRUN}
}


# ============modify_bashrc============
# Description: allow bashrc modification
# Usage: [[ $(modify_bashrc) == true ]]
#        modify_bashrc reset true
function modify_bashrc()
{
    if [[ $(dryrun) == true ]]; then
        echo false
        return
    fi    
    
    if [[ -z ${MODBASHRC} || "$1" == "reset" ]]; then
        if [[ -z ${2} ]]; then
            MODBASHRC=true
        else
            MODBASHRC=${2}
        fi
        return
    fi

    echo ${MODBASHRC}
}


# ============mexec============
# Description: Execute with root priviliged fro system wide installation
# Usage: mexec command
function mexec()
{
    if [[ $(install_systemwide) == true ]]; then
        sudo -E $@
    else
        $@
    fi
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
            #Remove with "dpkg -r yourpackagename"
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


# ============mapt-get============
# Description: Apt-get without prompt (ignores system wide setting)
# Usage: mapt-get install ...
function mapt-get()
{
    local pkgmgr=""
    if [[ $(cmdexists "apt-get") == true ]]; then
        pkgmgr="apt-get -y --allow-unauthenticated"
    elif [[ $(cmdexists "dnf") == true ]]; then
        pkgmgr="dnf -y"
    elif [[ $(cmdexists "yum") == true ]]; then
        pkgmgr="apt-get -y"
    elif [[ $(cmdexists "pkg") == true ]]; then
        pkgmgr="pkg -y"
    else
        pkgmgr="apt-get -y --allow-unauthenticated"
    fi

    require_web_access
    if [[ $(system_privileges) == true ]]; then
        sudo -E ${pkgmgr} "$@"
    else
        echo "Skip ${pkgmgr} $@ (no system priviliges)"
    fi
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


# ============addFile============
# Description: Add a line to a resource file when it doesn't exist yet
# Usage: addFile line filename
function addFile()
{
	local FILE="${1}"
	local LINE="${2}"

	grep "${LINE}" "${FILE}" > /dev/null 2>&1
	local Retval=$?
	if [[ $Retval != 0 ]]; then
	    cprint "Adding \"${LINE}\" to file ${FILE}"
        echo ${LINE} | mexec tee --append ${FILE}
        grep "${LINE}" "${FILE}" > /dev/null 2>&1
	    Retval=$?
	else
	    cprint "Line \"${LINE}\" already in ${FILE}"
	fi

	return $Retval
}


# ============addProfile============
# Description: Add a line to your environment file (e.g. .bashrc)
# Usage: addProfile filename line
function addProfile()
{
	addFile ${1} "${2}"

    if [[ $? == 0 && $(modify_bashrc) == true ]]; then
        if [[ $(install_systemwide) == true ]]; then
            addFile "/etc/bash.bashrc" "[ -r ${1} ] && source \"${1}\""
        else
            addFile "${HOME}/.bashrc" "[ -r ${1} ] && source \"${1}\""
        fi
    fi

    local Retval=$?
	return $Retval
}


# ============addVar============
# Description: Add variable
function addVar()
{
	export ${1}=${2}
}


# ============pathExists============
# Description: checks whether path exists
function pathExists()
{
    # Expand variables, commonly ${HOME}
    local _path=$(eval echo ${1})
    if [[ -e ${_path} ]];then
        echo true
    else
        echo false
    fi
}


# ============addBinPath============
# Description: Add path to ${PATH}
function addBinPath()
{
    if [[ $(pathExists ${1}) == false ]];then
        cprint "Does not exist: ${1}"
        return
    fi
	if [[ ":${PATH}:" != *":${1}:"* ]]; then
        export PATH=${1}:${PATH}
    fi
}


# ============addInclPath============
# Description: Add path to ${CPATH}
function addInclPath()
{
    if [[ $(pathExists ${1}) == false ]];then
        cprint "Does not exist: ${1}"
        return
    fi
	if [[ ":${CPATH}:" != *":${1}:"* ]]; then
        export CPATH=${1}:${CPATH}
    fi
}


# ============addLibPath============
# Description: Add path to ${LD_LIBRARY_PATH}
function addLibPath()
{
    if [[ $(pathExists ${1}) == false ]];then
        cprint "Does not exist: ${1}"
        return
    fi
    # For dynamic linking:
	if [[ ":${LD_LIBRARY_PATH},:" != *":${1}:"* ]]; then
        export LD_LIBRARY_PATH=${1}:${LD_LIBRARY_PATH}
    fi
    # For static linking:
	if [[ ":${LIBRARY_PATH},:" != *":${1}:"* ]]; then
        export LIBRARY_PATH=${1}:${LIBRARY_PATH}
    fi
}


# ============addPkgConfigPath============
# Description: Add path to ${PKG_CONFIG_PATH}
function addPkgConfigPath()
{
    if [[ $(pathExists ${1}) == false ]];then
        cprint "Does not exist: ${1}"
        return
    fi
	if [[ ":${PKG_CONFIG_PATH},:" != *":${1}:"* ]]; then
        export PKG_CONFIG_PATH=${1}:${PKG_CONFIG_PATH}
    fi
}


# ============addVarProfile============
# Description: Add variable premanently
function addVarProfile()
{
    addProfile ${1} "export ${2}=${3}"
}


# ============addBinPathProfile============
# Description: Add path to ${PATH} premanently
function addBinPathProfile()
{
    if [[ $(pathExists ${2}) == false ]];then
        cprint "Does not exist: ${2}"
        return
    fi
    addProfile ${1} "export PATH=${2}:\${PATH}"
}


# ============addInclPathProfile============
# Description: Add path to ${CPATH} premanently
function addInclPathProfile()
{
    if [[ $(pathExists ${2}) == false ]];then
        cprint "Does not exist: ${2}"
        return
    fi
    addProfile ${1} "export CPATH=${2}:\${CPATH}"
}


# ============addLibPathProfile============
# Description: Add path to ${LD_LIBRARY_PATH} premanently
function addLibPathProfile()
{
    if [[ $(pathExists ${2}) == false ]];then
        cprint "Does not exist: ${2}"
        return
    fi
    addProfile ${1} "export LD_LIBRARY_PATH=${2}:\${LD_LIBRARY_PATH}"
    addProfile ${1} "export LIBRARY_PATH=${2}:\${LIBRARY_PATH}"
}


# ============addPkgConfigPathProfile============
# Description: Add path to ${PKG_CONFIG_PATH} premanently
function addPkgConfigPathProfile()
{
    if [[ $(pathExists ${2}) == false ]];then
        cprint "Does not exist: ${2}"
        return
    fi
    addProfile ${1} "export PKG_CONFIG_PATH=${2}:\${PKG_CONFIG_PATH}"
}


# ============project_folder============
# Description: Project folder
function project_folder()
{
    echo "$( cd "$( dirname "${BASH_SOURCE[0]}" )"/../.. && pwd )"
}


# ============project_name============
# Description: Project name
function project_name()
{
    basename $(project_folder)
}


# ============project_echoresourcedir============
# Description:
function project_echoresourcedir()
{
    if [[ -z ${PROJECT_RESOURCE_DIR} ]]; then
        if [[ $(install_systemwide) == true ]]; then
            echo "/etc"
        else
            echo "${HOME}/"
        fi
    else
        echo ${PROJECT_RESOURCE_DIR} 
    fi
}


# ============project_resource============
# Description: Project resource file
function project_resource()
{
    if [[ $(install_systemwide) == true ]]; then
        echo "$(project_echoresourcedir)/$(project_name).bashrc"
    else
        echo "$(project_echoresourcedir)/.$(project_name)rc"
    fi
}


# ============project_echoprefix============
# Description:
function project_echoprefix()
{
    if [[ -z ${PROJECT_PREFIX} ]]; then
        echo ${1}
    else
        echo ${PROJECT_PREFIX}
    fi
}


# ============project_userbasestr============
# Description:
function project_userbasestr()
{
    project_echoprefix '${HOME}/.local'
}


# ============project_prefixstr============
# Description: 
function project_prefixstr()
{
    if [[ $(install_systemwide) == true ]]; then
        project_echoprefix "/usr/local"
    else
        project_userbasestr
    fi
}


# ============project_optstr============
# Description: 
function project_optstr()
{
    if [[ $(install_systemwide) == true ]]; then
        project_echoprefix "/opt"
    else
        project_prefixstr
    fi
}


# ============project_prefix============
# Description: use as configuration prefix so that files will be installed in
#                   $(project_prefix)/bin
#                   $(project_prefix)/lib
#                   ...
function project_prefix()
{
    fullpath $(project_prefixstr)
}


# ============project_opt============
# Description: use $(project_opt)/progname as configuration prefix so that files will be installed in
#                   $(project_opt)/progname/bin
#                   $(project_opt)/progname/lib
#                   ...
#              Don't forget to add $(project_opt)/progname/... to the environment:
#                   prefix=$(project_opt)/progname
#                   prefixstr=$(project_optstr)/progname
#                   addProfile $(project_resource) "# Installed progname: ${prefixstr}"
#                   addBinPath ${prefix}
#                   addBinPathProfile $(project_resource) ${prefixstr}
#                   addLibPath ${prefix}
#                   addLibPathProfile $(project_resource) ${prefixstr}
function project_opt()
{
    fullpath $(project_optstr)
}


# ============project_userbase============
# Description: see project_opt
function project_userbase()
{
    fullpath $(project_userbasestr)
}


# ============timer============
# Description: start and stop timer
function timer()
{
    if [[ -z ${START_TIME} || "$1" == "reset" ]]; then
        START_TIME=${SECONDS}
        return
    fi

    local elapsedtot=$((${SECONDS} - ${START_TIME}))
    local elapsedmin=$(( ${elapsedtot}/60 ))
    local elapsedsec=$(( ${elapsedtot}-60*${elapsedmin} ))
    cprint "Total execution time = ${elapsedmin} min ${elapsedsec} sec"
}


# ============cmdexists============
# Description: 
function cmdexists()
{
    if [[ -z $(command -v ${1}) ]]; then
        echo false
    else
        echo true
    fi
}


# ============libexists============
# Description: 
function libexists()
{
    for _path in ${LD_LIBRARY_PATH//:/ }; do
        if [[ -e ${_path}/${1}.so ]]; then
            echo true
            return
        fi
    done
    if [[ ! -z $(/sbin/ldconfig -p | grep -o -E "${1}\.so") ]]; then
        echo true
        return
    fi
    echo false
}


# ============libpath============
# Description: 
function libpath()
{
    for _path in ${LD_LIBRARY_PATH//:/ }; do
        if [[ -e ${_path}/${1}.so ]]; then
            echo ${_path}/${1}.so
            return
        fi
    done
    for _path in $(/sbin/ldconfig -p | grep -o -E "=> .+${1}\.so" | cut -c4- ); do
        if [[ -e ${_path} ]]; then
            echo ${_path}
            return
        fi
    done
}


# ============libversion============
# Description: 
function libversion()
{
    local _path=$(libpath ${1})
    if [[ -e ${_path} ]]; then
        readelf -d ${_path} | grep SONAME | grep -o -E "${1}\.so[\.0-9]+[0-9]" | grep -o -E "[0-9][\.0-9]*[0-9]?"
        return
    fi
    echo 0
}


# ============cmd_path============
# Description: 
function cmd_path()
{
    echo $(which ${1} | grep -E -o ".+/")
}


# ============cmd_full_bin============
# Description: 
function cmd_full_bin()
{
    echo $(type ${1}  | awk -v N=$3 '{print $3}')
}


# ============os_arch============
# Description: 
function os_arch()
{
	local osarch=$(getconf LONG_BIT)
	if [ $? -ne 0 ]; then
		osarch=$(uname -m | sed 's/x86_//;s/i[3-6]86/32/')
	fi
    echo ${osarch}
}


# ============install_info============
# Description: 
function install_info()
{
    cprint "Root priviliges: $(system_privileges)"
    cprint "System wide installation: $(install_systemwide)"
    cprint "Prefix for dependencies: $(project_prefix)"
    cprint "Opt directory: $(project_opt)"
    cprint "Resource file: $(project_resource)"
}

