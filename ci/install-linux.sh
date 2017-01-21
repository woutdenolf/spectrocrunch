#!/bin/bash

hcol='\033[0;35m'
ncol='\033[0m'

echo -e "${hcol}Download cmake ...${ncol}"
cd $CACHED_FOLDER
mkdir cmake
cd cmake
wget --no-check-certificate -q http://www.cmake.org/files/v3.7/cmake-3.7.2.tar.gz
tar xf cmake-3.7.2.tar.gz
cd cmake-3.7.2
echo -e "${hcol}Configure cmake ...${ncol}"
./configure
echo -e "${hcol}Build cmake ...${ncol}"
make -s
echo -e "${hcol}Install cmake ...${ncol}"
sudo make install -s

cd $TRAVIS_BUILD_DIR

echo -e "${hcol}Download SimpleElastix ...${ncol}"
cd $CACHED_FOLDER
mkdir simpleelastix
cd simpleelastix
git clone -b master https://github.com/kaspermarstal/SimpleElastix
mkdir build
cd build
echo -e "${hcol}Configure SimpleElastix ...${ncol}"
cmake ../SimpleElastix/SuperBuild 1>/dev/null
echo -e "${hcol}Build SimpleElastix ...${ncol}"
#travis_wait 30 mvn make 1>/dev/null
make -s
echo -e "${hcol}Install SimpleElastix ...${ncol}"
ls -R ./SimpleITK-build
cd ./SimpleITK-build/Wrapping/Python/Packaging
#cd ./SimpleITK-build/Wrapping/PythonPackage
python setup.py install

cd $TRAVIS_BUILD_DIR

