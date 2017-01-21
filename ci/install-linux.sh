#!/bin/bash

# Minimizing output:
#   1. longtask ...
#   2. make -s | grep -v '%]' | grep -v '^--' | ...
#   3. configure | grep -v '^--'

function longtask() {
    # $@: all parameters as a string
    # $!: PID of last job running in background.
    local cmd="$@"
    $cmd & 
    while true; do
        ps -p$! 1>/dev/null
        if [ $? = 0 ]; then
         echo "still running: $cmd"; sleep 1m
        else
         echo "$cmd finished"
         break
        fi
    done
}

hcol='\033[0;35m'
ncol='\033[0m'

#============cmake============

echo -e "${hcol}Download cmake ...${ncol}"
cd $CACHED_FOLDER
mkdir cmake
cd cmake
wget --no-check-certificate -q http://www.cmake.org/files/v3.7/cmake-3.7.2.tar.gz
tar -xzf cmake-3.7.2.tar.gz
cd cmake-3.7.2

echo -e "${hcol}Configure cmake ...${ncol}"
./configure | grep -v '^--'

echo -e "${hcol}Build cmake ...${ncol}"
make -s | grep -v '%]' | grep -v '^--'

echo -e "${hcol}Install cmake ...${ncol}"
sudo make install -s | grep -v '%]' | grep -v '^--'

#============SimpleElastix============

echo -e "${hcol}Download SimpleElastix ...${ncol}"
cd $CACHED_FOLDER
mkdir simpleelastix
cd simpleelastix
git clone -b master https://github.com/kaspermarstal/SimpleElastix
mkdir build
cd build

echo -e "${hcol}Configure SimpleElastix ...${ncol}"
longtask cmake -DCMAKE_RULE_MESSAGES=OFF -DCMAKE_INSTALL_MESSAGE=NEVER ../SimpleElastix/SuperBuild

echo -e "${hcol}Build SimpleElastix ...${ncol}"
longtask make -s | grep -v '%]' | grep -v '^--' | grep -v '^Installing'

echo -e "${hcol}Install SimpleElastix ...${ncol}"
ls -R ./SimpleITK-build
cd ./SimpleITK-build/Wrapping/Python/Packaging
#cd ./SimpleITK-build/Wrapping/PythonPackage
python setup.py install

#============cleanup============

cd $TRAVIS_BUILD_DIR

