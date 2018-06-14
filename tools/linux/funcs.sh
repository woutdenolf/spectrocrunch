#!/bin/bash
# 
# Helper functions.
# 

SCRIPT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $SCRIPT_ROOT/funcs-string.sh
source $SCRIPT_ROOT/funcs-essentials.sh

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
# Usage: [[ $(install_systemwide) ]]
#        install_systemwide reset true
function install_systemwide()
{
    if [[ -z ${INSTALL_SYSTEMWIDE} || "$1" == "reset" ]]; then
        if [[ -z ${2} ]]; then
            INSTALL_SYSTEMWIDE=$(system_privileges)
        else
            INSTALL_SYSTEMWIDE=${2}
        fi
    fi

    echo ${INSTALL_SYSTEMWIDE}
}


# ============dryrun============
# Description: check whether we are doing a dry run
# Usage: [[ $(dryrun) ]]
#        dryrun reset true
function dryrun()
{
    if [[ -z ${DRYRUN} || "$1" == "reset" ]]; then
        if [[ -z ${2} ]]; then
            DRYRUN=true
        else
            DRYRUN=${2}
        fi
    fi

    echo ${DRYRUN}
}


# ============mexec============
# Description: Execute with sudo if priviliged user
# Usage: mexec command
function mexec()
{
    if [[ $(system_privileges) == true ]]; then
        sudo -E $@
    else
        eval $@
    fi
}


# ============mmakeinstall============
# Description: Execute make install
# Usage: mmakeinstall
function mmakeinstall()
{
    if [[ $(system_privileges) == true && $(install_systemwide) == true ]]; then
        sudo -E checkinstall -y
    else
        make install -s
    fi
}


# ============mapt-get============
# Description: Apt-get without prompt
# Usage: mapt-get command
function mapt-get()
{
    require_web_access
    if [[ $(system_privileges) == true ]]; then
        sudo -E apt-get -y --force-yes $@
    else
        echo "Skip apt-get $@ (no system priviliges)"
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

    dpkg --compare-versions ${currentv} "le" ${requiredv}
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


# ============addFile============
# Description: Add a line to a resource file when it doesn't exist yet
# Usage: addFile line filename
function addFile()
{
	local FILE="${1}"
	local LINE="${2}"
	grep "${LINE}" "${FILE}" > /dev/null 2>&1
	Retval=$?
	if [ $Retval -ne 0 ]; then
		echo "${LINE}" >> "${FILE}"
        grep "${LINE}" "${FILE}" > /dev/null 2>&1
	    Retval=$?
	fi

	return $Retval
}


# ============addProfile============
# Description: Add a line to your environment file (e.g. .bashrc)
# Usage: addProfile line filename
function addProfile()
{
	addFile "$@"

    if [[ $? != 0 ]]; then
        if [[ $(install_systemwide) ]]; then
            addFile "/etc/bash.bashrc" "[ -r ${1} ] && source \"${1}\""
        else
            addFile "$HOME/.bashrc" "[ -r ${1} ] && source \"${1}\""
        fi
    fi

    Retval=$?
	return $Retval
}


# ============addBinPath============
# Description: Add path to $PATH
function addBinPath()
{
	if [[ ":$PATH:" != *":${1}:"* ]]; then
        export PATH=${1}:$PATH
    fi
}


# ============addInclPath============
# Description: Add path to $CPATH
function addInclPath()
{
	if [[ ":$CPATH:" != *":${1}:"* ]]; then
        export CPATH=${1}:$CPATH
    fi
}


# ============addLibPath============
# Description: Add path to $LD_LIBRARY_PATH
function addLibPath()
{
	if [[ ":$LD_LIBRARY_PATH,:" != *":${1}:"* ]]; then
        export LD_LIBRARY_PATH=${1}:$LD_LIBRARY_PATH
    fi
	if [[ ":$LIBRARY_PATH,:" != *":${1}:"* ]]; then
        export LIBRARY_PATH=${1}:$LIBRARY_PATH
    fi
}


# ============addBinPathProfile============
# Description: Add path to $PATH premanently
function addBinPathProfile()
{
    addProfile ${1} "export PATH=${2}:\$PATH"
}


# ============addInclPathProfile============
# Description: Add path to $CPATH premanently
function addInclPathProfile()
{
    addProfile ${1} "export CPATH=${2}:\$CPATH"
}


# ============addLibPathProfile============
# Description: Add path to $LD_LIBRARY_PATH premanently
function addLibPathProfile()
{
    addProfile ${1} "export LD_LIBRARY_PATH=${2}:\$LD_LIBRARY_PATH"
    addProfile ${1} "export LIBRARY_PATH=${2}:\$LIBRARY_PATH"
}


# ============project_resource============
# Description: Project name
function project_name()
{
    if [[ -z ${PROJECT} ]]; then
        PROJECT="project"
    fi

    echo ${PROJECT}
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


# ============timer============
# Description: 
function timer()
{
    if [[ -z ${START_TIME} || "$1" == "reset" ]]; then
        START_TIME=$SECONDS
    fi

    local ELAPSED_TIME=$((${SECONDS} - ${START_TIME}))
    cprint "Total execution time = $(( $ELAPSED_TIME/60 )) min"
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

