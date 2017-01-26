#!/bin/bash
# 
# This script will install all spectrocrunch Python 2 and 3 dependencies.
# 

# ============Check whether the current user is a sudoer============
if [[ ! -z "$((sudo -n true) 2>&1)" ]]; then
  echo $(whoami) "does not have sudo rights (login as root or any other sudoer)"
  return 1
fi

# ============Parse arguments============
show_help()
{
  echo "
        Usage: prepare_installation  [-v version] [-y] [-t] [-d]

        -v version      Python version to be used (2 or 3).
        -y              Answer yes to everything.
        -t              Time limited build.
        -d              Dry run

        For Example: ./prepare_installation -v 3 -d

        -h              Help
       "
}


unset CHOICE
OPTIND=0
PYTHONBIN="python"
PIPV="pip"
TIMELIMITED=false
TIMELEFT=true
NOTDRY=true
CI=false
while getopts "v:ythd" opt; do
  case $opt in
    h)
      show_help
      return 1
      ;;
    y)
      CHOICE="y"
      ;;
    t)
      TIMELIMITED=true
      ;;
    d)
      NOTDRY=false
      ;;
    v)
      if [ "$OPTARG" = "3" ]; then
        PYTHONBIN="python3"
        PIPV=$(compgen -c "pip-3" | sort | awk '{print $NF}')
      fi
      ;;
    \?)
      echo "Invalid option: -$OPTARG. Use -h flag for help." >&2
      return 1
      ;;
  esac
done
if [[ "$TRAVIS" == "true" ]]; then
    CHOICE="y"
    TIMELIMITED=true
    TIMELEFT=true
    CI=true
fi

# ============Initialize environment============
hcol='\033[0;36m'
ncol='\033[0m'

# Working directory for building dependencies: $(pwd)/2.7, $(pwd)/3.4, ...
RESTORE_WD=$(pwd)
PYTHONV=`$PYTHONBIN -c "import sys;t='{v[0]}.{v[1]}'.format(v=list(sys.version_info[:2]));sys.stdout.write(t)";`
mkdir -p $PYTHONV
cd $PYTHONV
INSTALL_WD=$(pwd)

echo -e "${hcol}Python version: $PYTHONV${ncol}"

PYTHON_EXECUTABLE=$(which $PYTHONBIN) # full path
PYTHON_INCLUDE_DIR="/usr/include/python$PYTHONV"
PYTHON_LIBRARY="/usr/lib/libpython$PYTHONV.so"

SPECTROCRUNCH_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"/..

if [[ "$(dnsdomainname)" == "esrf.fr" ]]; then
  echo -e "${hcol}Setting esrf proxy ...${ncol}"
  export http_proxy="http://proxy.esrf.fr:3128"
  export https_proxy="http://proxy.esrf.fr:3128"
fi

# ============Notifications============
if [[ -z ${CHOICE} ]]; then
  read -p "Approximately 12GB of data will added to \"$(pwd)\". Continue (Y/n)?" CHOICE
fi
case "$CHOICE" in 
  y|Y ) ;;
  n|N ) return 1;;
  * ) ;;
esac

# ============Install dependencies============
echo -e "${hcol}Installing basic ...${ncol}"
if [[ $NOTDRY == true ]]; then
  sudo -E apt-get -y install curl wget git
fi

echo -e "${hcol}Installing basic python ...${ncol}"
if [[ $NOTDRY == true ]]; then
  if [[ $CI == true ]]; then
    sudo -E apt-get -y install build-essential cmake $PYTHONBIN $PYTHONBIN-dev $PYTHONBIN-pip
  else
    sudo -E apt-get -y install build-essential cmake
  fi
fi

echo -e "${hcol}Installing python module dependencies ...${ncol}"
if [[ $NOTDRY == true ]]; then
  sudo -E apt-get -y install hdf5-devel # h5py
  sudo -E apt-get -y install libgeos-dev # shapely
  sudo -E apt-get -y install swig # xraylib
  sudo -E apt-get -y install opencl-headers $PYTHONBIN-pyopencl # silx
fi

# ============Install modules============
echo -e "${hcol}Installing python modules ...${ncol}"
if [[ $NOTDRY == true ]]; then
  sudo -E apt-get -y install $PYTHONBIN-numpy $PYTHONBIN-scipy $PYTHONBIN-matplotlib $PYTHONBIN-h5py $PYTHONBIN-setuptools $PYTHONBIN-pyopencl
  if [[ $CI == true ]]; then
    sudo -E -H $PIPV install --upgrade pip # problem on with Python 2 and 3 distros
    sudo -E -H $PIPV install --upgrade setuptools
  fi
  
  sudo -E -H $PIPV install -r $SPECTROCRUNCH_ROOT/requirements.txt
  ##sudo -E -H $PIPV install --upgrade ConfigParser # pyfdmnes (should be there by default)
fi

# ============Install xraylib============
if [ ! -f xraylib/xraylib-3.2.0/python/.libs/_xraylib.so ]; then
  echo -e "${hcol}Download xraylib ...${ncol}"
  mkdir -p xraylib
  cd xraylib
  if [[ $NOTDRY == true ]]; then
    curl -O http://lvserver.ugent.be/xraylib/xraylib-3.2.0.tar.gz
    tar -xvf xraylib-3.2.0.tar.gz
    cd xraylib-3.2.0
  fi

  echo -e "${hcol}Configure xraylib ...${ncol}"
  if [[ $NOTDRY == true ]]; then
    ./configure --enable-python --enable-python-integration PYTHON=$PYTHON_EXECUTABLE PYTHON_VERSION=$PYTHONV
  fi

  echo -e "${hcol}Build xraylib ...${ncol}"
  if [[ $NOTDRY == true ]]; then
    make -s -j2
    make check
  fi
else
  cd xraylib/xraylib-3.2.0
fi

echo -e "${hcol}Install xraylib ...${ncol}"
if [[ $NOTDRY == true ]]; then
  make install -s
fi

# ============Install fdmnes============
cd $INSTALL_WD
mkdir -p fdmnes
cd fdmnes

echo -e "${hcol}Install fdmnes ...${ncol}"
if [[ $NOTDRY == true ]]; then
  if [ -d /sware/exp/fdmnes ]; then
    if [ ! -d /opt/fdmnes ]; then
      ln -s /sware/exp/fdmnes /opt/fdmnes
    fi
    if [ ! -d /usr/local/bin/fdmnes ]; then
      ln -s /opt/fdmnes /usr/local/bin/fdmnes
    fi
  else
    curl -O http://neel.cnrs.fr/IMG/zip/fdmnes_2017_01_10.zip
    unzip fdmnes_2017_01_10.zip -d /usr/local
    ln -s /usr/local/fdmnes/fdmnes_linux64 /usr/local/bin/fdmnes
  fi
fi

echo -e "${hcol}Download pyfdmnes ...${ncol}"
if [[ $NOTDRY == true ]]; then
  git clone https://github.com/woutdenolf/pyFDMNES.git
fi

echo -e "${hcol}Install pyfdmnes ...${ncol}"
if [[ $NOTDRY == true ]]; then
  cd pyFDMNES
  $PYTHONBIN setup.py install
fi

# ============Install cmake > 3.0.0 (for spectrocrunch)============
dpkg --compare-versions "$(cmake --version | head -1 | awk '{print $3}')" "lt" "3.0.0"
if [ $? = 0 ]; then
  echo -e "${hcol}Download cmake ...${ncol}"
  cd $INSTALL_WD
  mkdir -p cmake
  cd cmake
  if [[ $NOTDRY == true ]]; then
    wget --no-check-certificate -q http://www.cmake.org/files/v3.7/cmake-3.7.2.tar.gz
    tar -xzf cmake-3.7.2.tar.gz
    cd cmake-3.7.2
  fi

  if [[ $TIMELEFT == true ]]; then
    echo -e "${hcol}Configure cmake ...${ncol}"
    if [[ $NOTDRY == true ]]; then
      ./configure
    fi

    echo -e "${hcol}Build cmake ...${ncol}"
    if [[ $NOTDRY == true ]]; then
      make -s -j2
    fi

    echo -e "${hcol}Install cmake ...${ncol}"
    if [[ $NOTDRY == true ]]; then
      sudo make install -s

      if [[ $TIMELIMITED == true ]]; then
        TIMELEFT=false
      fi
    fi
  fi
fi

# ============Install simpleelastix============
cd $INSTALL_WD
if [ ! -f simpleelastix/build/SimpleITK-build/Wrapping/Python/Packaging/setup.py ]; then
  echo -e "${hcol}Download SimpleElastix ...${ncol}"
  mkdir -p simpleelastix
  cd simpleelastix
  if [[ $NOTDRY == true ]]; then
    git clone https://github.com/kaspermarstal/SimpleElastix
  fi
  mkdir -p build
  cd build

  if [[ $TIMELEFT == true && ! -f Makefile ]]; then
    echo -e "${hcol}Configure SimpleElastix ...${ncol}"
    if [[ $NOTDRY == true ]]; then
      cmake ../SimpleElastix/SuperBuild -DPYTHON_EXECUTABLE:FILEPATH=$PYTHON_EXECUTABLE -DPYTHON_INCLUDE_DIR:PATH=$PYTHON_INCLUDE_DIR -DPYTHON_LIBRARY:FILEPATH=$PYTHON_LIBRARY
      if [[ $TIMELIMITED == true ]]; then
          TIMELEFT=false
      fi
    fi
  fi

  if [[ $TIMELEFT == true ]]; then
    echo -e "${hcol}Build SimpleElastix ...${ncol}"
    OMP_NUM_THREADS=2
    if [[ $NOTDRY == true ]]; then
      make -s -j2
    fi
  fi
else
  cd simpleelastix/build
fi

if [[ $TIMELEFT == true ]]; then
    echo -e "${hcol}Install SimpleElastix ...${ncol}"
    if [[ $NOTDRY == true ]]; then
      cd ./SimpleITK-build/Wrapping/Python/Packaging
      $PYTHONBIN setup.py install
    fi
fi

# ============Cleanup============
echo -e "${hcol}Cleaning up ...${ncol}"
cd $RESTORE_WD

if [[ $NOTDRY == true ]]; then
  apt-get -y autoremove

  if [[ $TIMELEFT == true ]]; then
      echo -e "${hcol}All done! You should not be able to install spectrocrunch.${ncol}"
  else
      echo -e "${hcol}Not everything has been build due to time restrictions. Run the script again.${ncol}"
  fi
fi

