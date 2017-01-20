#!/bin/bash

cwd=$(pwd)
echo $cwd

echo "Installing cmake ..."
wget --no-check-certificate http://www.cmake.org/files/v3.2/cmake-3.2.2.tar.gz
tar xf cmake-3.2.2.tar.gz
cd cmake-3.2.2
./configure
make
make install
cd $cwd

echo "Installing SimpleElastix ..."
mkdir simpleelastix
cd simpleelastix
git clone -b master https://github.com/kaspermarstal/SimpleElastix
mkdir build
cd build
cmake ../SimpleElastix/SuperBuild
make
cd ./SimpleITK-build/Wrapping/Python/Packaging
python setup.py install
cd $cwd

