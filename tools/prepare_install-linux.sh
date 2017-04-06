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
        -d              Dry run

        For Example: ./prepare_installation -v 3 -d

        -h              Help
       "
}

# ============Initialize environment============
START_TIME=$SECONDS

hcol='\033[0;36m'
ncol='\033[0m'

RESTORE_WD=$(pwd)
export SPECTROCRUNCH_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"/..
if [[ -z "$((sudo -n true) 2>&1)" ]]; then
  export SYSTEM_PRIVILIGES=true 
else
  export SYSTEM_PRIVILIGES=false
fi

RETURNCODE_ARG=1
RETURNCODE_PYTHONENV=2
RETURNCODE_CANCEL=3
export BUILDSTEPS=0
export BUILDSTEP=0

if [[ "$(dnsdomainname)" == "esrf.fr" ]]; then
  echo -e "${hcol}Setting esrf proxy ...${ncol}"
  export http_proxy="http://proxy.esrf.fr:3128"
  export https_proxy="http://proxy.esrf.fr:3128"
fi

# ============Parse arguments============
export PYTHONBINAPT="python"
export PIPBINAPT="pip"
export TIMELIMITED=false
export TIMELEFT=true
export NOTDRY=true

OPTIND=0
while getopts "v:ythd" opt; do
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

# ============Python============
echo -e "${hcol}Looking for python interpreter ...${ncol}"

export PYTHONBIN=$PYTHONBINAPT

if [[ -z `which $PYTHONBIN` && $SYSTEM_PRIVILIGES == true && $NOTDRY == true ]]; then
  sudo -E apt-get install $PYTHONBINAPT $PYTHONBINAPT-dev $PYTHONBINAPT-qt4
fi

if [ -z `which $PYTHONBIN` ]; then
  echo -e "${hcol}$PYTHONBIN is not installed on this system.${ncol}"
  cd $RESTORE_WD
  return $RETURNCODE_PYTHONENV
fi

# ============Check python version============
PYTHONMAJORV=`$PYTHONBIN -c "import sys;print(sys.version_info[0])";`
PYTHONV=`$PYTHONBIN -c "import sys;t='{v[0]}.{v[1]}'.format(v=list(sys.version_info[:2]));print(t)";`

if [[ $PYTHONMAJORV == "3" ]]; then
  dpkg --compare-versions "${PYTHONV}" "lt" "3.4"
  if [ $? = 0 ]; then
    echo -e "${hcol}Python version must be >= 3.4 (used ${PYTHONV}).${ncol}"
    cd $RESTORE_WD
    return $RETURNCODE_PYTHONENV
  fi
else
  dpkg --compare-versions "${PYTHONV}" "lt" "2.7"
  if [ $? = 0 ]; then
    echo -e "${hcol}Python version must be >= 2.7 (used ${PYTHONV}).${ncol}"
    cd $RESTORE_WD
    return $RETURNCODE_PYTHONENV
  fi
fi

mkdir -p ${PYTHONV}
cd ${PYTHONV}
INSTALL_WD=$(pwd)

# ============Pip============
echo -e "${hcol}Looking for pip package ...${ncol}"

export PIPBIN=$PIPBINAPT

if [[ -z `which $PIPBIN` && $SYSTEM_PRIVILIGES == true && $NOTDRY == true ]]; then
  sudo -E apt-get install $PIPBINAPT
fi

if [ -z `which $PIPBIN` ]; then
  echo -e "${hcol}$PIPBIN is not installed on this system.${ncol}"
  cd $RESTORE_WD
  return $RETURNCODE_PYTHONENV
fi

echo -e "${hcol}Upgrading pip ...${ncol}"
if [[ $NOTDRY == true ]]; then
  $PIPBIN install --upgrade pip
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

PYTHON_EXECUTABLE=$(which $PYTHONBIN) # full path
PYTHONFULLV=`$PYTHONBIN -c "import sys;t='{v[0]}.{v[1]}.{v[2]}'.format(v=list(sys.version_info[:3]));print(t)";`
PYTHON_INCLUDE_DIR=`$PYTHONBIN -c "import distutils.sysconfig; print(distutils.sysconfig.get_python_inc());"`
PYTHON_LIBRARY=`$PYTHONBIN -c "import distutils.sysconfig,os; print(os.path.join(distutils.sysconfig.get_config_var('LIBDIR'),distutils.sysconfig.get_config_var('LDLIBRARY')));"`

echo -e "${hcol}Python version: $PYTHONFULLV ${ncol}"
echo -e "${hcol}Python location: $PYTHON_EXECUTABLE ${ncol}"
echo -e "${hcol}Python include: $PYTHON_INCLUDE_DIR ${ncol}"
echo -e "${hcol}Python library: $PYTHON_LIBRARY ${ncol}"
echo -e "${hcol}Pip:$($PIPBIN --version| awk '{$1= ""; print $0}')${ncol}"

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
echo -e "${hcol}Install basics ...${ncol}"
if [[ $NOTDRY == true && $SYSTEM_PRIVILIGES == true ]]; then
  sudo -E apt-get -y install build-essential cmake curl wget git
fi

BUILDSTEP=$(( $BUILDSTEP+1 ))
BUILDSTEPS=$(( $BUILDSTEPS+1 ))

# ============Install system dependencies============
echo -e "${hcol}Install python module dependencies ...${ncol}"
if [[ $NOTDRY == true && $SYSTEM_PRIVILIGES == true ]]; then
  sudo -E apt-get install $PYTHONBINAPT-qt4 # pymca
  sudo -E apt-get -y install libgeos-dev # shapely
  sudo -E apt-get -y install opencl-headers # pyopencl
  sudo -E apt-get -y install libffi-dev # pyopencl
  sudo -E apt-get -y install libgl1-mesa-dev libglu1-mesa-dev mesa-common-dev # pymca
fi

BUILDSTEP=$(( $BUILDSTEP+1 ))
BUILDSTEPS=$(( $BUILDSTEPS+1 ))

# ============Install modules============
echo -e "${hcol}Install python modules ...${ncol}"
if [[ $NOTDRY == true ]]; then
  $PIPBIN install --upgrade setuptools
  $PIPBIN install --upgrade wheel
  $PIPBIN install --upgrade numpy # silx
  $PIPBIN install --upgrade mako # pyopencl

  $PIPBIN install --upgrade -r $SPECTROCRUNCH_ROOT/requirements.txt
  $PIPBIN install --upgrade --egg pymca #TODO: wait for pymca to get fixed
fi

BUILDSTEP=$(( $BUILDSTEP+1 ))
BUILDSTEPS=$(( $BUILDSTEPS+1 ))

# ============Custom installation============
source $SPECTROCRUNCH_ROOT/tools/install-xraylib.sh
cd $INSTALL_WD

source $SPECTROCRUNCH_ROOT/tools/install-fdmnes.sh
cd $INSTALL_WD

source $SPECTROCRUNCH_ROOT/tools/install-simpleelastix.sh
cd $INSTALL_WD

# ============Cleanup============
echo -e "${hcol}Cleaning up ...${ncol}"
cd $RESTORE_WD

if [[ $NOTDRY == true ]]; then
  if [[ $SYSTEM_PRIVILIGES == true ]]; then
    sudo -E apt-get -y autoremove
  else
    PATH=$HOME/.local:$PATH
    echo -e "${hcol}Make sure to add PATH=\$HOME/.local:\$PATH to your bashrc file.${ncol}"
  fi

  if [[ $TIMELEFT == true ]]; then
      echo -e "${hcol}All done ($BUILDSTEP/$BUILDSTEPS)! You should now be able to install spectrocrunch.${ncol}"
  else
      echo -e "${hcol}Not everything has been build due to time restrictions. Run the script again ($BUILDSTEP/$BUILDSTEPS).${ncol}"
  fi
else
  echo -e "${hcol}Dry build $BUILDSTEP/$BUILDSTEPS.${ncol}"
fi

ELAPSED_TIME=$(($SECONDS - $START_TIME))
echo -e "${hcol}Total execution time = $(( $ELAPSED_TIME/60 )) min${ncol}"


