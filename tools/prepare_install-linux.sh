#!/bin/bash
# 
# This script will install all spectrocrunch Python 2 and 3 dependencies.
# 

# ============Usage============
show_help()
{
  echo "
        Usage: prepare_installation  [-v version] [-y] [-t] [-d]

        -v version      Python version to be used (2, 3, 2.7, 3.5, ...).
        -y              Answer yes to everything.
        -t              Time limited build.
        -d              Dry run.
        -u              Install for user only.

        For Example: ./prepare_installation -v 3 -d

        -h              Help
       "
}

# ============Initialize environment============
SCRIPT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $SCRIPT_ROOT/funcs.sh
resetEnv

# ============Adapt environment based on script arguments============
OPTIND=0
while getopts "v:uythd" opt; do
  case $opt in
    h)
      show_help
      cd $RESTORE_WD
      return $RETURNCODE_ARG
      ;;
    y)
      FORCECHOICE="y"
      ;;
    t)
      TIMELIMITED=true
      ;;
    u)
      SYSTEM_PRIVILIGES=false
      ;;
    d)
      NOTDRY=false
      ;;
    \?)
      echo "Invalid option: -$OPTARG. Use -h flag for help." >&2
      cd $RESTORE_WD
      return $RETURNCODE_ARG
      ;;
  esac
done
initEnv

# ============Initialize Python============
initPython
Retval=$?
if [ $Retval -ne 0 ]; then
    cd $RESTORE_WD
    return $Retval
fi

mkdir -p ${PYTHONV}
cd ${PYTHONV}
INSTALL_WD=$(pwd)

# ============Initialize Pip============
initPip
Retval=$?
if [ $Retval -ne 0 ]; then
    cd $RESTORE_WD
    return $Retval
fi

# ============Show information and ask to proceed============
# /usr/lib/python2.7/dist-packages: system-wide, installed with package manager
# /usr/local/lib/python2.7/dist-packages: system-wide, installed with pip
#
# /usr/lib/python2.7/site-packages:
# /usr/local/lib/python2.7/site-packages:
#
# /usr/lib/python2.7:
# /usr/lib/pymodules/python2.7:
#
# ~/.local/lib/python2.7/site-packages: user-only, installed with pip --user
#
# <virtualenv_name>/lib/python2.7/site-packages: virtual environment, installed with pip

cprint "Python version: $PYTHONFULLV"
cprint "Python location: $PYTHON_EXECUTABLE"
cprint "Python include: $PYTHON_INCLUDE_DIR"
cprint "Python library: $PYTHON_LIBRARY"
cprint "Pip:$($PIPBIN --version| awk '{$1= ""; print $0}')"

if [[ -z $FORCECHOICE ]]; then
    read -p "Approximately 12GB of data will added to \"$(pwd)\". Continue (Y/n)?" CHOICE
else
    CHOICE=$FORCECHOICE
fi
case "$CHOICE" in 
  y|Y ) ;;
  n|N ) 
        cd $RESTORE_WD
        return $RETURNCODE_CANCEL;;
  * ) ;;
esac

# ============Install basics============
cprint "Install basics ..."
if [[ $NOTDRY == true && $SYSTEM_PRIVILIGES == true ]]; then
    sudo -E apt-get -y install make build-essential cmake curl wget git
fi

BUILDSTEP=$(( $BUILDSTEP+1 ))
BUILDSTEPS=$(( $BUILDSTEPS+1 ))

# ============Install system dependencies============
cprint "Install python module dependencies ..."
if [[ $SYSTEM_PRIVILIGES == true ]]; then
    if [[ $NOTDRY == true ]]; then
        sudo -E apt-get -y install $PYTHONBINAPT-qt4 # pymca
        sudo -E apt-get -y install libgeos-dev # shapely
        sudo -E apt-get -y install opencl-headers # pyopencl
        sudo -E apt-get -y install libffi-dev # pyopencl
        sudo -E apt-get -y install libgl1-mesa-dev libglu1-mesa-dev mesa-common-dev # pymca
    fi
    BUILDSTEP=$(( $BUILDSTEP+1 ))
    BUILDSTEPS=$(( $BUILDSTEPS+1 ))
else
    source $SCRIPT_ROOT/install-opencl.sh # pyopencl
    cd $INSTALL_WD

    source $SCRIPT_ROOT/install-libgeos.sh # shapely
    cd $INSTALL_WD
fi

# ============Install modules============
cprint "Install python modules ..."
if [[ $NOTDRY == true ]]; then
    $PIPBIN install --upgrade setuptools
    $PIPBIN install --upgrade wheel
    $PIPBIN install --upgrade numpy # silx
    $PIPBIN install --upgrade mako # pyopencl

    $PIPBIN install --upgrade -r $SCRIPT_ROOT/../requirements.txt
    $PIPBIN install --upgrade --egg pymca #TODO: wait for pymca to get fixed
fi

BUILDSTEP=$(( $BUILDSTEP+1 ))
BUILDSTEPS=$(( $BUILDSTEPS+1 ))

# ============Custom installation============
source $SCRIPT_ROOT/install-xraylib.sh
cd $INSTALL_WD

source $SCRIPT_ROOT/install-fdmnes.sh
cd $INSTALL_WD

source $SCRIPT_ROOT/install-simpleelastix.sh
cd $INSTALL_WD

# ============Cleanup============
cprint "Cleaning up ..."
cd $RESTORE_WD

if [[ $NOTDRY == true ]]; then
    if [[ $SYSTEM_PRIVILIGES == true ]]; then
        sudo -E apt-get -y autoremove
    else
        cprint "Variables have been added to $SPECTROCRUNCHRC."
    fi

    if [[ $TIMELEFT == true ]]; then
        cprint "All done ($BUILDSTEP/$BUILDSTEPS)! You should now be able to install spectrocrunch."
    else
        cprint "Not everything has been build due to time restrictions. Run the script again ($BUILDSTEP/$BUILDSTEPS)."
    fi
else
    cprint "Dry build $BUILDSTEP/$BUILDSTEPS."
fi

ELAPSED_TIME=$(($SECONDS - $START_TIME))
cprint "Total execution time = $(( $ELAPSED_TIME/60 )) min"


