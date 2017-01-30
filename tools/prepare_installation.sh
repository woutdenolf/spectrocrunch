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
PYTHONBINAPT="python"
PIPBINAPT="pip"
PYTHONBIN=$PYTHONBINAPT
PIPBIN=$PIPBINAPT
PYTHONMAJORV=2
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
      if [ "$(echo $OPTARG | head -c 1)" = "3" ]; then
        PYTHONMAJORV=3
        PYTHONBINAPT="python3"
        PIPBINAPT="pip3"
        PYTHONBIN="python$OPTARG"
        PIPBIN="pip$OPTARG"
      elif [ "$(echo $OPTARG | head -c 1)" = "2" ]; then
        PYTHONMAJORV=2
        PYTHONBINAPT="python"
        PIPBINAPT="pip"
        PYTHONBIN="python$OPTARG"
        PIPBIN="pip$OPTARG"
      else
        echo "Python version number must start with 2 or 3." >&2
        return 1
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

RESTORE_WD=$(pwd)
SPECTROCRUNCH_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"/..

if [[ "$(dnsdomainname)" == "esrf.fr" ]]; then
  echo -e "${hcol}Setting esrf proxy ...${ncol}"
  export http_proxy="http://proxy.esrf.fr:3128"
  export https_proxy="http://proxy.esrf.fr:3128"
fi

# ============Install basics============
echo -e "${hcol}Installing basics ...${ncol}"
if [[ $NOTDRY == true ]]; then
  sudo -E apt-get -y install curl wget git
  sudo -E apt-get -y install build-essential cmake $PYTHONBINAPT $PYTHONBINAPT-dev $PYTHONBINAPT-pip
fi

# ============Python version============
function setpythonbin() {
  if [[ "${PYTHONBIN#*.}" != "$PYTHONBIN" ]]; then
    # Already a specific version like 2.7
    return
  fi

  PYTHONBIN=$(compgen -c "python$1" | grep "^python$1\..$" | sort | tail -n -1)
  if [[ -z $PYTHONBIN ]]; then
    PYTHONBIN=$(compgen -c "python-$1" | grep "^python-$1\..$" | sort | tail -n -1)
  fi
  if [[ -z $PYTHONBIN ]]; then
    PYTHONBIN=$PYTHONBINAPT
  fi
}
setpythonbin $PYTHONMAJORV
PYTHON_EXECUTABLE=$(which $PYTHONBIN) # full path

PYTHONV=`$PYTHONBIN -c "import sys;t='{v[0]}.{v[1]}'.format(v=list(sys.version_info[:2]));sys.stdout.write(t)";`
PYTHONFULLV=`$PYTHONBIN -c "import sys;t='{v[0]}.{v[1]}.{v[2]}'.format(v=list(sys.version_info[:3]));sys.stdout.write(t)";`

PYTHON_INCLUDE_DIR=`$PYTHONBIN -c "import distutils.sysconfig; print(distutils.sysconfig.get_python_inc());"`

PYTHON_LIBRARY_NAME=`$PYTHONBIN -c "import distutils.sysconfig; print(distutils.sysconfig.get_config_var('LDLIBRARY'));"`
PYTHON_PREFIX=="/usr/"
PYTHON_LIBRARY="${PYTHON_PREFIX}lib/${PYTHON_LIBRARY_NAME}"
if [[ ! -f $PYTHON_LIBRARY ]]; then
  PYTHON_PREFIX="/usr/local/"
  PYTHON_LIBRARY="${PYTHON_PREFIX}lib/${PYTHON_LIBRARY_NAME}"
fi

if [[ $PYTHONMAJORV = 3 ]]; then
  dpkg --compare-versions "${PYTHONV}" "lt" "3.4"
  if [ $? = 0 ]; then
    echo -e "${hcol}Python version must be >= 3.4 (used ${PYTHONV}).${ncol}"
    return 1
  fi
else
  dpkg --compare-versions "${PYTHONV}" "lt" "2.7"
  if [ $? = 0 ]; then
    echo -e "${hcol}Python version must be >= 2.7 (used ${PYTHONV}).${ncol}"
    return 1
  fi
fi

mkdir -p ${PYTHONV}
cd ${PYTHONV}
INSTALL_WD=$(pwd)

# ============Pip version============
function setpipbin() {
  if [[ "${PIPBIN#*.}" != "$PIPBIN" ]]; then
    # Already a specific version like 2.7
    return
  fi

  PIPBIN=$(compgen -c "pip$1" | grep "^pip$1\..$" | sort | tail -n -1)
  if [[ -z $PIPBIN ]]; then
    PIPBIN=$(compgen -c "pip-$1" | grep "^pip-$1\..$" | sort | tail -n -1)
  fi
  if [[ -z $PIPBIN ]]; then
    PIPBIN=$PIPBINAPT
  fi
}
setpipbin $PYTHONMAJORV
echo -e "${hcol}Upgrading pip ...${ncol}"
sudo -E -H $PIPBIN install --upgrade pip

# ============Notifications============
echo -e "${hcol}Python version: $PYTHONFULLV ($PYTHON_EXECUTABLE)${ncol}"
echo -e "${hcol}Pip:$($PIPBIN --version| awk '{$1= ""; print $0}')${ncol}"

if [[ -z ${CHOICE} ]]; then
  read -p "Approximately 12GB of data will added to \"$(pwd)\". Continue (Y/n)?" CHOICE
fi
case "$CHOICE" in 
  y|Y ) ;;
  n|N ) return 1;;
  * ) ;;
esac

# ============Install dependencies============
echo -e "${hcol}Installing python module dependencies ...${ncol}"
if [[ $NOTDRY == true ]]; then
  sudo -E apt-get -y install hdf5-devel # h5py
  sudo -E apt-get -y install libgeos-dev # shapely
  sudo -E apt-get -y install swig # xraylib, simpleelastix
  sudo -E apt-get -y install opencl-headers # pyopencl
  sudo -E apt-get -y install libffi-dev # pyopencl
  sudo -E apt-get -y install libinsighttoolkit-dev # simpleelastix
fi

# ============Install modules============
echo -e "${hcol}Installing python modules ...${ncol}"
if [[ $NOTDRY == true ]]; then
  
  # Try to install modules with apt-get first
  sudo -E apt-get -y install $PYTHONBINAPT-numpy
  sudo -E apt-get -y install $PYTHONBINAPT-scipy
  sudo -E apt-get -y install $PYTHONBINAPT-h5py
  sudo -E apt-get -y install $PYTHONBINAPT-setuptools
  sudo -E apt-get -y install $PYTHONBINAPT-pyopencl
  sudo -E apt-get -y install $PYTHONBINAPT-matplotlib
  sudo -E apt-get -y install $PYTHONBINAPT-pyopencl

  
  sudo -E -H $PIPBIN install --upgrade setuptools
  sudo -E -H $PIPBIN install --upgrade numpy # silx
  sudo -E -H $PIPBIN install --upgrade mako # pyopencl

  setpipbin $PYTHONMAJORV

  sudo -E -H $PIPBIN install --upgrade -r $SPECTROCRUNCH_ROOT/requirements.txt
  sudo -E -H $PIPBIN install --upgrade pymca #TODO: doesn't want to build
fi

# ============Install xraylib============
if [ ! -f xraylib/xraylib-3.2.0/python/.libs/_xraylib.so ]; then
  echo -e "${hcol}Download xraylib ...${ncol}"
  mkdir -p xraylib
  cd xraylib
  if [[ $NOTDRY == true && ! -d xraylib-3.2.0 ]]; then
    curl -O http://lvserver.ugent.be/xraylib/xraylib-3.2.0.tar.gz
    tar -xvf xraylib-3.2.0.tar.gz
    cd xraylib-3.2.0
  fi

  echo -e "${hcol}Configure xraylib ...${ncol}"
  if [[ $NOTDRY == true ]]; then
    ./configure --enable-python --enable-python-integration PYTHON=$PYTHON_EXECUTABLE PYTHON_VERSION=$PYTHONV PYTHON_PREFIX=$PYTHON_PREFIX
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
  sudo -E make install -s
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
    if [[ ! -f /usr/local/bin/fdmnes ]]; then
      curl -O http://neel.cnrs.fr/IMG/zip/fdmnes_2017_01_10.zip
      sudo -E unzip fdmnes_2017_01_10.zip -d /usr/local
      sudo -E ln -s /usr/local/fdmnes/fdmnes_linux64 /usr/local/bin/fdmnes
    fi
  fi
fi

echo -e "${hcol}Download pyfdmnes ...${ncol}"
if [[ $NOTDRY == true && ! -d pyFDMNES ]]; then
  git clone https://github.com/woutdenolf/pyFDMNES.git
fi

echo -e "${hcol}Install pyfdmnes ...${ncol}"
if [[ $NOTDRY == true ]]; then
  cd pyFDMNES
  $PYTHONBIN setup.py build -f
  $PYTHONBIN setup.py install
fi

# ============Install cmake > 3.0.0 (for spectrocrunch)============
dpkg --compare-versions "$(cmake --version | head -1 | awk '{print $3}')" "lt" "3.0.0"
if [ $? = 0 ]; then
  echo -e "${hcol}Download cmake ...${ncol}"
  cd $INSTALL_WD
  mkdir -p cmake
  cd cmake
  if [[ $NOTDRY == true && ! -d cmake-3.7.2 ]]; then
    wget --no-check-certificate http://www.cmake.org/files/v3.7/cmake-3.7.2.tar.gz
    tar -xzf cmake-3.7.2.tar.gz
  fi
  cd cmake-3.7.2

  if [[ $TIMELEFT == true && ! -f Makefile ]]; then
    echo -e "${hcol}Configure cmake ...${ncol}"
    if [[ $NOTDRY == true ]]; then
      ./configure
    fi

    echo -e "${hcol}Build cmake ...${ncol}"
    if [[ $NOTDRY == true ]]; then
      make -s -j2
    fi

    if [[ $TIMELIMITED == true ]]; then
      TIMELEFT=false
    fi
  fi

  echo -e "${hcol}Install cmake ...${ncol}"
  if [[ $NOTDRY == true ]]; then
    sudo make install -s
  fi
fi

# ============Install simpleelastix============
cd $INSTALL_WD
if [ ! -f simpleelastix/build/SimpleITK-build/Wrapping/Python/Packaging/setup.py ]; then
  echo -e "${hcol}Download SimpleElastix ...${ncol}"
  mkdir -p simpleelastix
  cd simpleelastix
  if [[ $NOTDRY == true && ! -d SimpleElastix ]]; then
    git clone https://github.com/kaspermarstal/SimpleElastix
  fi
  mkdir -p build
  cd build

  if [[ $TIMELEFT == true && ! -f Makefile ]]; then
    echo -e "${hcol}Configure SimpleElastix ...${ncol}"
    if [[ $NOTDRY == true ]]; then
      # TODO: run twice to get the right python interpreter?
      # TODO: USE_SYSTEM_ITK, USE_SYSTEM_ELASTIX
      cmake -DBUILD_EXAMPLES:BOOL=OFF \
            -DBUILD_SHARED_LIBS:BOOL=OFF \
            -DBUILD_TESTING:BOOL=OFF \
            -DUSE_SYSTEM_SWIG:BOOL=ON \
            -DPYTHON_EXECUTABLE:FILEPATH=$PYTHON_EXECUTABLE \
            -DPYTHON_INCLUDE_DIR:PATH=$PYTHON_INCLUDE_DIR \
            -DPYTHON_LIBRARY:FILEPATH=$PYTHON_LIBRARY \
            -DWRAP_CSHARP:BOOL=OFF \
            -DWRAP_JAVA:BOOL=OFF \
            -DWRAP_LUA:BOOL=OFF \
            -DWRAP_PYTHON:BOOL=ON \
            -DWRAP_R:BOOL=OFF \
            -DWRAP_RUBY:BOOL=OFF \
            -DWRAP_TCL:BOOL=OFF \
            ../SimpleElastix/SuperBuild
      cmake -DBUILD_EXAMPLES:BOOL=OFF \
            -DBUILD_SHARED_LIBS:BOOL=OFF \
            -DBUILD_TESTING:BOOL=OFF \
            -DUSE_SYSTEM_SWIG:BOOL=ON \
            -DPYTHON_EXECUTABLE:FILEPATH=$PYTHON_EXECUTABLE \
            -DPYTHON_INCLUDE_DIR:PATH=$PYTHON_INCLUDE_DIR \
            -DPYTHON_LIBRARY:FILEPATH=$PYTHON_LIBRARY \
            -DWRAP_CSHARP:BOOL=OFF \
            -DWRAP_JAVA:BOOL=OFF \
            -DWRAP_LUA:BOOL=OFF \
            -DWRAP_PYTHON:BOOL=ON \
            -DWRAP_R:BOOL=OFF \
            -DWRAP_RUBY:BOOL=OFF \
            -DWRAP_TCL:BOOL=OFF \
            ../SimpleElastix/SuperBuild

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
      if [[ $TIMELIMITED == true ]]; then
          TIMELEFT=false
      fi
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
  sudo -E apt-get -y autoremove

  if [[ $TIMELEFT == true ]]; then
      echo -e "${hcol}All done! You should not be able to install spectrocrunch.${ncol}"
  else
      echo -e "${hcol}Not everything has been build due to time restrictions. Run the script again.${ncol}"
  fi
fi

