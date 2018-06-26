#!/bin/bash
# 
# This script will install pre-builds on Travis.
# 

# Download pre-build libraries
PYTHONV=`python -c "import sys;t='{v[0]}.{v[1]}.{v[2]}'.format(v=list(sys.version_info[:3]));print(t)"`

cd ${BUILD_FOLDER}

if [ ! -d ${PYTHONV} ]; then
    FILE=spectrocrunch.travis.python${PYTHONV}.tgz
    LINK1=https://transfer.sh/12avMO/${FILE}
    LINK2=http://ftp.esrf.fr/tmp/${FILE}
    MARKER=${BUILD_FOLDER}/${FILE}.done

    # Download to cache folder
    cd ${CACHED_FOLDER}
    ls -all
    if [ ! -f ${MARKER} ]; then
        echo "Download pre-build libraries ..."
        wget ${LINK1}
        ls -all
        if [ ! -f ${FILE} ]; then
            wget ${LINK2}
        fi

        # Unpack in build folder
        if [ -f ${FILE} ]; then
            echo "Unpack pre-build libraries ..."
            tar -xzf ${FILE} -C ${BUILD_FOLDER}
            touch ${MARKER}
            rm -f ${FILE}
        fi
    fi

    cd ${BUILD_FOLDER}
fi

# List pre-build libraries   
if [ -d ${PYTHONV} ]; then
    echo "Pre-build libraries:"
    ls ${PYTHONV}
else
    echo "No pre-build libraries"
fi

# Display when needed
export DISPLAY=:99.0
sh -e /etc/init.d/xvfb start

sudo -E add-apt-repository universe

