#!/bin/bash

#============definitions============

function reduceoutput() {
    # Usage:
    #   reduceoutput make
    #
    # Reducing output can also be done with ... | grep -v '%]' | grep -v '^--'

    # Travis needs output every 10 minutes, otherwise it will shut down
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
        echo -e "${hcol}Status $1 is already ${!1}${ncol}"
    else
        export "$1"=0
        echo -e "${hcol}Status $1 is initialized to ${ncol}"
    fi 
}

function incstate() {
    let "$1++"
    echo -e "${hcol}Status $1 is incremented to ${!1}${ncol}"
}

function incandsavestate() {
    incstate $1
    savestate $1
}

function restoresavedstate() {
    if [ -f $ENV_FILE ]; then
        source $ENV_FILE
    fi    
}

function clearsavedstate() {
    if [ -f $ENV_FILE ]; then
        rm $ENV_FILE
    fi    
}

#============init============

restoresavedstate
initstate STATE_CMAKE
initstate STATE_SIMPLEELASTIX

#============cmake============

cd $CACHED_FOLDER
STATE=1

if [ "$STATE_CMAKE" -lt "$STATE" ]; then
    echo -e "${hcol}Download cmake ...${ncol}"
    mkdir -p cmake
    cd cmake
    wget --no-check-certificate -q http://www.cmake.org/files/v3.7/cmake-3.7.2.tar.gz
    tar -xzf cmake-3.7.2.tar.gz
    cd cmake-3.7.2

    echo -e "${hcol}Configure cmake ...${ncol}"
    ./configure

    incandsavestate STATE_CMAKE
else
    cd cmake/cmake-3.7.2
fi

incstate STATE

if [ "$STATE_CMAKE" -lt "$STATE" ]; then
    echo -e "${hcol}Build cmake ...${ncol}"
    make -s

    incandsavestate STATE_CMAKE
fi

incstate STATE

if [ "$STATE_CMAKE" -lt "$STATE" ]; then
    echo -e "${hcol}Install cmake ...${ncol}"
    sudo make install -s

    incandsavestate STATE_CMAKE
fi

#============SimpleElastix============

cd $CACHED_FOLDER
STATE=1

if [ "$STATE_SIMPLEELASTIX" -lt "$STATE" ]; then
    echo -e "${hcol}Download SimpleElastix ...${ncol}"
    mkdir -p simpleelastix
    cd simpleelastix
    git clone -b master https://github.com/kaspermarstal/SimpleElastix
    mkdir -p build
    cd build

    incandsavestate STATE_SIMPLEELASTIX
else
    cd simpleelastix/build
fi

incstate STATE

OMP_NUM_THREADS=2
if [ "$STATE_SIMPLEELASTIX" -lt "$STATE" ]; then
    echo -e "${hcol}Configure SimpleElastix ...${ncol}"
    cmake -DCMAKE_RULE_MESSAGES=OFF -DCMAKE_INSTALL_MESSAGE=NEVER -DBUILD_TESTING=OFF ../SimpleElastix/SuperBuild

    incandsavestate STATE_SIMPLEELASTIX
fi

incstate STATE

if [ "$STATE_SIMPLEELASTIX" -lt "$STATE" ]; then
    echo -e "${hcol}Build SimpleElastix ...${ncol}"
    make -s -j2

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
clearsavedstate
savestate STATE_CMAKE
savestate STATE_SIMPLEELASTIX

cd $TRAVIS_BUILD_DIR

