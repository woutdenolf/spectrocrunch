#!/bin/bash

echo "Installing SimpleElastix ..."

cwd=$(pwd)
echo $cwd
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

