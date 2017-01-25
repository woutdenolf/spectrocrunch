#!/bin/bash
# 
# This script will install all spectrocrunch dependencies
# 
#    prepare_installation [-y]
#

# ============Check whether the current user is a sudoer============
if [[ ! -z "$((sudo -n true) 2>&1)" ]]; then
  echo $(whoami) "does not have sudo rights (login as root or any other sudoer)"
  return 1
fi

# ============Parse arguments============
unset choice
OPTIND=0
while getopts ":y" opt; do
  case $opt in
    y)
      choice="y"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      ;;
  esac
done

# ============Notifications============
if [[ -z ${choice} ]]; then
  read -p "Approximately 11GB of data will be downloaded to $(pwd). Continue (Y/n)?" choice
fi
case "$choice" in 
  y|Y ) ;;
  n|N ) return 1;;
  * ) ;;
esac

# ============Initialize environment============
hcol='\033[0;36m'
ncol='\033[0m'
INSTALL_ROOT=$(pwd)
SPECTROCRUNCH_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"/..
if [[ "$(dnsdomainname)" == "esrf.fr" ]]; then
  echo -e "${hcol}Setting esrf proxy ...${ncol}"
  export http_proxy="http://proxy.esrf.fr:3128"
  export https_proxy="http://proxy.esrf.fr:3128"
fi

# ============Install dependencies============
echo -e "${hcol}Installing basic ...${ncol}"
sudo -E apt-get -y install curl wget git

echo -e "${hcol}Installing basic python ...${ncol}"
sudo -E apt-get -y install build-essential cmake python python-dev python-pip

echo -e "${hcol}Installing python module dependencies ...${ncol}"
sudo -E apt-get -y install hdf5-devel # h5py
sudo -E apt-get -y install libgeos-dev # shapely
sudo -E apt-get -y install swig # xraylib
sudo -E apt-get -y install opencl-headers python-pyopencl # silx

# ============Install modules============
sudo -E apt-get -y install python-numpy python-scipy python-matplotlib python-h5py python-setuptools python-pyopencl
sudo -E -H pip install -r $SPECTROCRUNCH_ROOT/requirements.txt

# ============Install xraylib============
if [ ! -f xraylib/xraylib-3.2.0/python/.libs/_xraylib.so ]; then
  echo -e "${hcol}Download xraylib ...${ncol}"
  mkdir xraylib
  cd xraylib
  curl -O http://lvserver.ugent.be/xraylib/xraylib-3.2.0.tar.gz
  tar -xvf xraylib-3.2.0.tar.gz
  cd xraylib-3.2.0

  echo -e "${hcol}Configure xraylib ...${ncol}"
  ./configure --enable-python --enable-python-integration

  echo -e "${hcol}Build xraylib ...${ncol}"
  make -s
  make check
else
  cd xraylib/xraylib-3.2.0
fi

echo -e "${hcol}Install xraylib ...${ncol}"
make install -s

cd $INSTALL_ROOT

# ============Install fdmnes============
mkdir fdmnes
cd fdmnes

echo -e "${hcol}Install fdmnes ...${ncol}"
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

echo -e "${hcol}Download pyfdmnes ...${ncol}"
git clone https://github.com/woutdenolf/pyFDMNES.git

echo -e "${hcol}Install pyfdmnes ...${ncol}"
cd pyFDMNES
python setup.py install

cd $INSTALL_ROOT

# ============Install cmake > 3.0.0 (for spectrocrunch)============
dpkg --compare-versions "$(cmake --version | head -1 | awk '{print $3}')" "lt" "3.0.0"
if [ $? = 0 ]; then
  echo -e "${hcol}Download cmake ...${ncol}"
  mkdir -p cmake
  cd cmake
  wget --no-check-certificate -q http://www.cmake.org/files/v3.7/cmake-3.7.2.tar.gz
  tar -xzf cmake-3.7.2.tar.gz
  cd cmake-3.7.2

  echo -e "${hcol}Configure cmake ...${ncol}"
  ./configure

  echo -e "${hcol}Build cmake ...${ncol}"
  make -s

  echo -e "${hcol}Install cmake ...${ncol}"
  sudo make install -s
  cd $INSTALL_ROOT
fi

# ============Install simpleelastix============
if [ ! -f simpleelastix/build/SimpleITK-build/Wrapping/Python/Packaging/setup.py ]; then
  echo -e "${hcol}Download SimpleElastix ...${ncol}"
  mkdir -p simpleelastix
  cd simpleelastix
  git clone https://github.com/kaspermarstal/SimpleElastix
  mkdir -p build
  cd build

  echo -e "${hcol}Configure SimpleElastix ...${ncol}"
  #cmake ../SimpleElastix/SuperBuild

  echo -e "${hcol}Build SimpleElastix ...${ncol}"
  OMP_NUM_THREADS=2
  #make -s -j2
else
 cd simpleelastix/build
fi

echo -e "${hcol}Install SimpleElastix ...${ncol}"
#cd ./SimpleITK-build/Wrapping/Python/Packaging
#python setup.py install

cd $INSTALL_ROOT

# ============Cleanup============
echo -e "${hcol}Cleaning up ...${ncol}"
apt-get -y autoremove
echo -e "${hcol}All done! You should not be able to install spectrocrunch.${ncol}"

