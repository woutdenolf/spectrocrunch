#!/bin/bash
# 
# Helper functions.
# 

# ============addFile============
# Add a line to a resource file
addFile()
{
	local FILE="$1"
	local LINE="$2"
	grep "$LINE" "$FILE" > /dev/null 2>&1
	Retval=$?
	if [ $Retval -ne 0 ]; then
		echo "$LINE" >> "$FILE"
        grep "$LINE" "$FILE" > /dev/null 2>&1
	    Retval=$?
	fi

	return $Retval
}


# ============addProfile============
# Add a line to your environment file (e.g. .bashrc)
addProfile()
{
	addFile "$@"

    if [[ $1 == $SPECTROCRUNCHRC && $SYSTEM_PRIVILIGES == false ]]; then
        addFile "$HOME/.bashrc" "[ -r $SPECTROCRUNCHRC ] && source \"$SPECTROCRUNCHRC\""
    fi

    Retval=$?
	return $Retval
}

# ============addBinPath============
# Add path to $PATH
addBinPath()
{
	if [[ ":$PATH:" != *":$1:"* ]]; then
        export PATH=$1:$PATH
    fi
}

# ============addInclPath============
# Add path to $CPATH
addInclPath()
{
	if [[ ":$CPATH:" != *":$1:"* ]]; then
        export CPATH=$1:$CPATH
    fi
}


# ============addLibPath============
# Add path to $LD_LIBRARY_PATH
addLibPath()
{
	if [[ ":$LD_LIBRARY_PATH,:" != *":$1:"* ]]; then
        export LD_LIBRARY_PATH=$1:$LD_LIBRARY_PATH
    fi
	if [[ ":$LIBRARY_PATH,:" != *":$1:"* ]]; then
        export LIBRARY_PATH=$1:$LIBRARY_PATH
    fi
}

# ============addBinPathProfile============
# Add path to $PATH premanently
addBinPathProfile()
{
    addProfile $1 "export PATH=$2:\$PATH"
}

# ============addInclPathProfile============
# Add path to $CPATH premanently
addInclPathProfile()
{
    addProfile $1 "export CPATH=$2:\$CPATH"
}

# ============addLibPathProfile============
# Add path to $LD_LIBRARY_PATH premanently
addLibPathProfile()
{
    addProfile $1 "export LD_LIBRARY_PATH=$2:\$LD_LIBRARY_PATH"
    addProfile $1 "export LIBRARY_PATH=$2:\$LIBRARY_PATH"
}

# ============initPython============
initPython()
{
    echo -e "${hcol}Looking for python interpreter ...${ncol}"

    export PYTHONBIN=$PYTHONBINAPT

    if [[ -z `which $PYTHONBIN` && $SYSTEM_PRIVILIGES == true && $NOTDRY == true ]]; then
        sudo -E apt-get -y install $PYTHONBINAPT $PYTHONBINAPT-dev $PYTHONBINAPT-qt4
    fi

    if [[ -z `which $PYTHONBIN` ]]; then
        echo -e "${hcol}$PYTHONBIN is not installed on this system.${ncol}"
        return $RETURNCODE_PYTHONENV
    fi

    initEnv

    # Check version
    if [[ $PYTHONMAJORV == "3" ]]; then
        dpkg --compare-versions "${PYTHONV}" "lt" "3.4"
        if [ $? = 0 ]; then
            echo -e "${hcol}Python version must be >= 3.4 (used ${PYTHONV}).${ncol}"
            return $RETURNCODE_PYTHONENV
        fi
    else
        dpkg --compare-versions "${PYTHONV}" "lt" "2.7"
        if [ $? = 0 ]; then
            echo -e "${hcol}Python version must be >= 2.7 (used ${PYTHONV}).${ncol}"
            return $RETURNCODE_PYTHONENV
        fi
    fi

    return $RETURNCODE_OK
}

# ============initPip============
initPip()
{
    echo -e "${hcol}Looking for pip package ...${ncol}"

    export PIPBIN=$PIPBINAPT

    if [[ -z `which $PIPBIN` && $SYSTEM_PRIVILIGES == true && $NOTDRY == true ]]; then
        sudo -E apt-get -y install $PIPBINAPT
    fi

    if [[ -z `which $PIPBIN` ]]; then
        echo -e "${hcol}$PIPBIN is not installed on this system.${ncol}"
        return $RETURNCODE_PYTHONENV
    fi

    echo -e "${hcol}Upgrading pip ...${ncol}"
    if [[ $NOTDRY == true ]]; then
        $PIPBIN install --upgrade pip
    fi

    initEnv

    return $RETURNCODE_OK
}

# ============_initEnv============
# Initialize common variables
# export: available to subprocesses
_initEnv()
{
    # Reset 
    local RESET=false
    if [ $# -ge 1 ]; then
        if [[ "$1" == true ]]; then
            RESET=true
        fi
    fi

    # ============Constants============
    RETURNCODE_OK=0
    RETURNCODE_ARG=1
    RETURNCODE_PYTHONENV=2
    RETURNCODE_CANCEL=3

    # ============Platform============
	OS_ARCH=$(getconf LONG_BIT)
	if [ $? -ne 0 ]; then
		OS_ARCH=$(uname -m | sed 's/x86_//;s/i[3-6]86/32/')
	fi
	OS_DISTRO=$(lsb_release -si)
	OS_VERSION=$(lsb_release -sr)

    # ============Python============
    if [[ -z $PYTHONBIN || $RESET == true ]]; then
        export PYTHONBIN=python
    fi
    if [[ -z $PYTHONBINAPT || $RESET == true ]]; then
        export PYTHONBINAPT=python
    fi
    if [[ -z $PIPBIN || $RESET == true ]]; then
        export PIPBIN=pip
    fi
    if [[ -z $PIPBINAPT || $RESET == true ]]; then
        export PIPBINAPT=pip
    fi

    PYTHONMAJORV=`$PYTHONBIN -c "import sys;print(sys.version_info[0])";`
    PYTHONV=`$PYTHONBIN -c "import sys;t='{v[0]}.{v[1]}'.format(v=list(sys.version_info[:2]));print(t)";`
    PYTHONFULLV=`$PYTHONBIN -c "import sys;t='{v[0]}.{v[1]}.{v[2]}'.format(v=list(sys.version_info[:3]));print(t)";`

    PYTHON_EXECUTABLE=$(which $PYTHONBIN) # full path
    
    PYTHON_INCLUDE_DIR=`$PYTHONBIN -c "import distutils.sysconfig; print(distutils.sysconfig.get_python_inc());"`
    PYTHON_LIBRARY=`$PYTHONBIN -c "import distutils.sysconfig,os; print(os.path.join(distutils.sysconfig.get_config_var('LIBDIR'),distutils.sysconfig.get_config_var('LDLIBRARY')));"`

    # ============User============
    if [[ -z $SYSTEM_PRIVILIGES || $RESET == true ]]; then
        if [[ -z "$((sudo -n true) 2>&1)" ]]; then
            SYSTEM_PRIVILIGES=true 
        else
            SYSTEM_PRIVILIGES=false
        fi
    fi

    if [[ -z $SPECTROCRUNCHRC || $RESET == true ]]; then
        if [[ $SYSTEM_PRIVILIGES == true ]]; then
            SPECTROCRUNCHRC="/etc/profile.d/spectrocrunchrc.sh"
        else
            SPECTROCRUNCHRC="$HOME/.spectrocrunchrc"
        fi
    fi

    if [[ -z $SPECTROCRUNCHLOCAL || $RESET == true ]]; then
        if [[ $SYSTEM_PRIVILIGES == true ]]; then
            SPECTROCRUNCHLOCAL="/usr/local"
        else
            SPECTROCRUNCHLOCAL="$HOME/.local"
        fi
    fi

    if [[ -z $SPECTROCRUNCHLOCALSTR || $RESET == true ]]; then
        if [[ $SYSTEM_PRIVILIGES == true ]]; then
            SPECTROCRUNCHLOCALSTR="/usr/local"
        else
            SPECTROCRUNCHLOCALSTR="\$HOME/.local"
        fi
    fi

    if [[ -z $SPECTROCRUNCHOPT || $RESET == true ]]; then
        if [[ $SYSTEM_PRIVILIGES == true ]]; then
            SPECTROCRUNCHOPT="/opt"
        else
            SPECTROCRUNCHOPT="$HOME/.local"
        fi
    fi

    if [[ -z $SPECTROCRUNCHOPTSTR || $RESET == true ]]; then
        if [[ $SYSTEM_PRIVILIGES == true ]]; then
            SPECTROCRUNCHOPTSTR="/opt"
        else
            SPECTROCRUNCHOPTSTR="\$HOME/.local"
        fi
    fi

    # ============Installation progress============
    if [[ -z $NOTDRY || $RESET == true ]]; then
        NOTDRY=true
    fi

    if [[ -z $BUILDSTEP || $RESET == true ]]; then
        BUILDSTEP=0
        BUILDSTEPS=0
    fi

    if [[ -z $TIMELEFT || $RESET == true ]]; then
        TIMELEFT=true
    fi

    if [[ -z $TIMELIMITED || $RESET == true ]]; then
        TIMELIMITED=false
    fi

    if [[ -z $START_TIME || $RESET == true ]]; then
        START_TIME=$SECONDS
    fi

    if [[ -z $RESTORE_WD || $RESET == true ]]; then
        RESTORE_WD=$(pwd)
    fi

    if [[ -z $INSTALL_WD || $RESET == true ]]; then
        INSTALL_WD=$(pwd)
    fi

    # ============Site specific============
    if [[ "$(dnsdomainname)" == "esrf.fr" ]]; then
        if [[ -z $http_proxy ]]; then
            export http_proxy="http://proxy.esrf.fr:3128"
        fi
        if [[ -z $http_proxy ]]; then
            export https_proxy="http://proxy.esrf.fr:3128"
        fi
    fi
}

# ============initEnv============
# Initialize common variables (only overwrite the derived variables)
initEnv()
{
    _initEnv false
}

# ============reinitEnv============
# Reset common variables (overwrite all)
resetEnv()
{
    _initEnv true
}

# ============reinitEnv============
# Output in color
cprint()
{
    local hcol='\033[0;36m'
    local ncol='\033[0m'
    echo -e "${hcol}$1 ${ncol}"
}

# ============mexec============
# Execute with sudo if priviliged user
mexec()
{
    if [[ $SYSTEM_PRIVILIGES == true ]]; then
        sudo -E $1
    else
        eval $1
    fi
}


