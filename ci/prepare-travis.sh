#!/bin/bash
# 
# This script will install all spectrocrunch Python 2 and 3 dependencies for CI.
# 

# Download pre-build libraries
PYTHONV=`python -c "import sys;t='{v[0]}.{v[1]}'.format(v=list(sys.version_info[:2]));print(t)";`

cd $BUILD_FOLDER

if [ ! -d ${PYTHONV} ]; then
    FILE=spectrocrunch.travis.python${PYTHONV}.tgz
    LINK1=https://transfer.sh/12avMO/$FILE
    LINK2=http://ftp.esrf.fr/tmp/$FILE

    # Download to cache folder
    cd $CACHED_FOLDER
    if [ ! -f $FILE ]; then
        echo "Download pre-build libraries ..."
        wget ${LINK1}
        if [ -f $FILE ]; then
            wget ${LINK2}
        fi
    fi

    # Unpack in build folder
    if [ -f $FILE ]; then
        echo "Unpack pre-build libraries ..."
        tar -xzf $FILE -C $BUILD_FOLDER
    fi
fi

# List pre-build libraries
cd $BUILD_FOLDER
ls -all
    
if [ -d ${PYTHONV} ]; then
    echo "Pre-build libraries:"
    ls ${PYTHONV}
else
    echo "No pre-build libraries"
fi

