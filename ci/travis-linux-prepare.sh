#!/bin/bash
# 
# This script will prepare Travis.
# 

function travis_download_prebuild()
{
    # Download pre-build libraries
    local PYTHONV=`python -c "import sys;t='{v[0]}.{v[1]}.{v[2]}'.format(v=list(sys.version_info[:3]));print(t)"`

    local DEP_FOLDER=dep_${PYTHONV}

    if [[ ! -d ${DEP_FOLDER} ]]; then
        local FILE=spectrocrunch.travis.python${PYTHONV}.tgz
        local LINK1=http://ftp.esrf.fr/tmp/${FILE}
        local LINK2=https://transfer.sh/12avMO/${FILE}
        
        # Download to cache folder
        if [[ ! -d ${DEP_FOLDER} ]]; then
            echo "Download pre-build libraries ..."
            wget ${LINK1}
            if [[ ! -f ${FILE} ]]; then
                wget ${LINK2}
            fi

            # Unpack in build folder
            if [[ -f ${FILE} ]]; then
                echo "Unpack pre-build libraries ..."
                tar -xzf ${FILE}
                rm -f ${FILE}
                sudo -E chown -R $(id -un):$(id -gn) ${DEP_FOLDER}
            fi
        fi
    fi

    # List pre-build libraries  
    if [[ -d ${DEP_FOLDER} ]]; then
        echo "Pre-build libraries:"
        pwd
        ls ${DEP_FOLDER}/*
    else
        echo "No pre-build libraries"
        pwd
        ls -all
    fi
}

function main()
{
    travis_download_prebuild

    # Display when needed
    export DISPLAY=:99.0
    sudo -E apt-get install xvfb
    sudo -E chmod +x /etc/init.d/xvfb
    sudo -E /etc/init.d/xvfb start

    # Add package repositories
    sudo -E add-apt-repository universe
    sudo -E apt-key update
    sudo -E apt-get update
}

main

