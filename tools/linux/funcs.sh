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
    echo -e "${hcol}${1} ${ncol}"
}


# ============cerror============
# Description: output to stdout in color
# Usage: cerror "..."
function cerror()
{
    local hcol='\033[0;31m'
    local ncol='\033[0m'
    echo -e "${hcol}${1} ${ncol}"
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
        return 
    fi

    if [[ $(system_privileges) == false ]]; then
        echo false
        return
    fi

    echo ${INSTALL_SYSTEMWIDE}
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
        return
    fi

    echo ${DRYRUN}
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
    if [[ $(install_systemwide) == true ]]; then
        local name=${1}
        if [[ -z ${name} ]];then
            name=$(randomstring 6)
        fi
        sudo -E checkinstall -y --pkgname "${name}-checkinstall"
    else
        make install -s
    fi
}


# ============mapt-get============
# Description: Apt-get without prompt (ignores system wide setting)
# Usage: mapt-get install ...
function mapt-get()
{
    local pkgmgr=""
    if [[ $(cmdexists "apt-get") == true ]]; then
        pkgmgr="apt-get -y --force-yes"
    elif [[ $(cmdexists "dnf") == true ]]; then
        pkgmgr="dnf -y"
    elif [[ $(cmdexists "yum") == true ]]; then
        pkgmgr="apt-get -y"
    elif [[ $(cmdexists "pkg") == true ]]; then
        pkgmgr="pkg -y"
    else
        pkgmgr="apt-get -y --force-yes"
    fi

    require_web_access
    if [[ $(system_privileges) == true ]]; then
        sudo -E ${pkgmgr} "$@"
    else
        echo "Skip ${pkgmgr} $@ (no system priviliges)"
    fi
}


# ============require_new_version============
# Description: we require a new version when current < required
# Usage: [[ $(require_new_version currentversion requiredversion) ]]
function require_new_version()
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

    dpkg --compare-versions ${currentv} "lt" ${requiredv}
    if [[ $? == 0 ]]; then
        echo true # current < required
    else
        echo false # current >= required
    fi
}


# ============require_new_version_strict============
# Description: 
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
        echo ${LINE} | mexec tee --append ${FILE}
        grep "${LINE}" "${FILE}" > /dev/null 2>&1
	    Retval=$?
	fi

	return $Retval
}


# ============addProfile============
# Description: Add a line to your environment file (e.g. .bashrc)
# Usage: addProfile filename line
function addProfile()
{
	addFile ${1} "${2}"

    if [[ $? == 0 ]]; then
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


# ============addBinPath============
# Description: Add path to ${PATH}
function addBinPath()
{
	if [[ ":${PATH}:" != *":${1}:"* ]]; then
        export PATH=${1}:${PATH}
    fi
}


# ============addInclPath============
# Description: Add path to ${CPATH}
function addInclPath()
{
	if [[ ":${CPATH}:" != *":${1}:"* ]]; then
        export CPATH=${1}:${CPATH}
    fi
}


# ============addLibPath============
# Description: Add path to ${LD_LIBRARY_PATH}
function addLibPath()
{
	if [[ ":${LD_LIBRARY_PATH},:" != *":${1}:"* ]]; then
        export LD_LIBRARY_PATH=${1}:${LD_LIBRARY_PATH}
    fi
	if [[ ":${LIBRARY_PATH},:" != *":${1}:"* ]]; then
        export LIBRARY_PATH=${1}:${LIBRARY_PATH}
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
    addProfile ${1} "export PATH=${2}:\${PATH}"
}


# ============addInclPathProfile============
# Description: Add path to ${CPATH} premanently
function addInclPathProfile()
{
    addProfile ${1} "export CPATH=${2}:\${CPATH}"
}


# ============addLibPathProfile============
# Description: Add path to ${LD_LIBRARY_PATH} premanently
function addLibPathProfile()
{
    addProfile ${1} "export LD_LIBRARY_PATH=${2}:\${LD_LIBRARY_PATH}"
    addProfile ${1} "export LIBRARY_PATH=${2}:\${LIBRARY_PATH}"
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


# ============project_resource============
# Description: Project resource file
function project_resource()
{
    if [[ $(install_systemwide) == true ]]; then
        echo "/etc/$(project_name).bashrc"
    else
        echo "${HOME}/.$(project_name)rc"
    fi
}


# ============project_prefixstr============
# Description: 
function project_prefixstr()
{
    if [[ $(install_systemwide) == true ]]; then
        echo "/usr/local"
    else
        echo "\${HOME}/.local"
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


# ============project_optstr============
# Description: 
function project_optstr()
{
    if [[ $(install_systemwide) == true ]]; then
        echo "/opt"
    else
        project_prefixstr
    fi
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


# ============project_userbasestr============
# Description:
function project_userbasestr()
{
    echo "\${HOME}/.local"
}


# ============project_userbase============
# Description: see project_opt
function project_userbase()
{
    fullpath $(project_userbasestr)
}


# ============timer============
# Description: 
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
    echo osarch
}


# ============install_info============
# Description: 
function install_info()
{
    cprint "Root priviliges: $(system_privileges)"
    cprint "System wide installation: $(install_systemwide)"
    cprint "Prefix for dependencies: $(project_prefix)"
    cprint "Opt directory: $(project_opt)"
}

