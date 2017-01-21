#!/bin/bash

cwd=$(pwd)
echo $cwd
echo $TRAVIS_BUILD_DIR
echo $HOME

hcol='\033[0;35m'
ncol='\033[0m'

echo -e "${hcol}Download cmake ...${ncol}"
cd $HOME
wget --no-check-certificate http://www.cmake.org/files/v3.7/cmake-3.7.2.tar.gz
tar xf cmake-3.7.2.tar.gz
cd cmake-3.7.2
echo -e "${hcol}Configure cmake ...${ncol}"
./configure 1>/dev/null
echo -e "${hcol}Build cmake ...${ncol}"
make 1>/dev/null
echo -e "${hcol}Install cmake ...${ncol}"
sudo make install 1>/dev/null

cd $cwd

echo -e "${hcol}Download SimpleElastix ...${ncol}"
cd $HOME
mkdir simpleelastix
cd simpleelastix
git clone -b master https://github.com/kaspermarstal/SimpleElastix
mkdir build
cd build
echo -e "${hcol}Configure SimpleElastix ...${ncol}"
cmake ../SimpleElastix/SuperBuild 1>/dev/null
echo -e "${hcol}Build SimpleElastix ...${ncol}"
travis_wait 30 mvn make 1>/dev/null
echo -e "${hcol}Install SimpleElastix ...${ncol}"
cd ./SimpleITK-build/Wrapping/Python/Packaging
python setup.py install 1>/dev/null

cd $cwd

