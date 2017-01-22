#!/bin/bash

#============definitions============

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

export ENV_FILE=$CACHED_FOLDER/env.dat
export hcol='\033[0;35m'
export ncol='\033[0m'

function savestate() {
    echo "export" "$1"="${!1}" >> $ENV_FILE
}

function initstate() {
    if [[ -v $1 ]]; then
        echo -e "${hcol}Status $1 is currently ${!1}${ncol}"
    else
        export "$1"=0
        echo -e "${hcol}Status $1 is initialized ${ncol}"
    fi 
}

function incstate() {
    let "$1++"
    
}

function incandsavestate() {
    let "$1++"
    savestate $1
}

function restorestate() {
    if [ -f $ENV_FILE ]; then
        source $ENV_FILE
    fi    
}

function clearstate() {
    if [ -f $ENV_FILE ]; then
        rm $ENV_FILE
    fi    
}

#============init============

restorestate
initstate STATE_CMAKE
initstate STATE_SIMPLEELASTIX

#============cmake============

cd $CACHED_FOLDER
STATE=1

if [ "$STATE_CMAKE" -lt "$STATE" ]; then
    echo -e "${hcol}Download cmake ...${ncol}"
    mkdir cmake
    cd cmake
    wget --no-check-certificate -q http://www.cmake.org/files/v3.7/cmake-3.7.2.tar.gz
    tar -xzf cmake-3.7.2.tar.gz

    cd cmake-3.7.2
    echo -e "${hcol}Configure cmake ...${ncol}"
    ./configure | grep -v '^--'

    incandsavestate STATE_CMAKE
else
    cd cmake/cmake-3.7.2
fi

incstate STATE

if [ "$STATE_CMAKE" -lt "$STATE" ]; then
    echo -e "${hcol}Build cmake ...${ncol}"
    make -s | grep -v '%]' | grep -v '^--'

    incandsavestate STATE_CMAKE
fi

incstate STATE

if [ "$STATE_CMAKE" -lt "$STATE" ]; then
    echo -e "${hcol}Install cmake ...${ncol}"
    sudo make install -s | grep -v '%]' | grep -v '^--'

    incandsavestate STATE_CMAKE
fi

#============SimpleElastix============

cd $CACHED_FOLDER
STATE=1

if [ "$STATE_SIMPLEELASTIX" -lt "$STATE" ]; then
    echo -e "${hcol}Download SimpleElastix ...${ncol}"
    cd $CACHED_FOLDER
    mkdir simpleelastix
    cd simpleelastix
    git clone -b master https://github.com/kaspermarstal/SimpleElastix
    mkdir build
    cd build

    incandsavestate STATE_SIMPLEELASTIX
else
    cd simpleelastix/build
fi

incstate STATE

OMP_NUM_THREADS=2
if [ "$STATE_SIMPLEELASTIX" -lt "$STATE" ]; then
    echo -e "${hcol}Configure SimpleElastix ...${ncol}"
    longtask cmake -DCMAKE_RULE_MESSAGES=OFF -DCMAKE_INSTALL_MESSAGE=NEVER -DBUILD_TESTING=OFF ../SimpleElastix/SuperBuild

    incandsavestate STATE_SIMPLEELASTIX
fi

incstate STATE

if [ "$STATE_SIMPLEELASTIX" -lt "$STATE" ]; then
    echo -e "${hcol}Build SimpleElastix ...${ncol}"
    #longtask make -s | grep -v '%]' | grep -v '^--' | grep -v '^Installing'
    make -s -j2 | grep -v '%]'

    incandsavestate STATE_SIMPLEELASTIX
fi

incstate STATE

if [ "$STATE_SIMPLEELASTIX" -lt "$STATE" ]; then
    echo -e "${hcol}Install SimpleElastix ...${ncol}"
    cd ./SimpleITK-build/Wrapping/Python/Packaging
    #cd ./SimpleITK-build/Wrapping/PythonPackage
    python setup.py install

    incandsavestate STATE_SIMPLEELASTIX
fi

#============cleanup============
clearstate
savestate STATE_CMAKE
savestate STATE_SIMPLEELASTIX

cd $TRAVIS_BUILD_DIR

