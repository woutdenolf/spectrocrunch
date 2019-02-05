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

function virtualdisplay()
{
    export DISPLAY=:99.0
    sudo -E apt-get install xvfb
    if [[ -f /etc/init.d/xvfb ]];then
        sudo -E chmod +x /etc/init.d/xvfb
        sudo -E /etc/init.d/xvfb start
    else
        local filename="/etc/systemd/system/xvfb.service"
        if [[ ! -f ${filename} ]];then
            sudo -E echo "[Unit]" > ${filename}
            sudo -E echo "Description=X Virtual Frame Buffer Service" >> ${filename}
            sudo -E echo "After=network.target" >> ${filename}
            sudo -E echo "" >> ${filename}
            sudo -E echo "[Service]" >> ${filename}
            sudo -E echo "ExecStart=/usr/bin/Xvfb ${DISPLAY} -screen 0 1024x768x24" >> ${filename}
            sudo -E echo "" >> ${filename}
            sudo -E echo "[Install]" >> ${filename}
            sudo -E echo "WantedBy=multi-user.target" >> ${filename}
        fi
        sudo -E systemctl enable ${filename}
        sudo -E service xvfb start
    fi
}

function main()
{
    travis_download_prebuild

    # Display
    virtualdisplay
    
    # Add package repositories
    #sudo -E add-apt-repository universe
    #sudo -E apt-key update
    #sudo -E apt-get update
}

main

