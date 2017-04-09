#!/bin/bash
# 
# This script will install all spectrocrunch Python 2 and 3 dependencies for CI.
# 

# Download pre-build libraries
PYTHONV=`python -c "import sys;t='{v[0]}.{v[1]}'.format(v=list(sys.version_info[:2]));print(t)";`

if [ ! -d ${PYTHONV} ]; then
    FILE=travis.python${PYTHONV}.tgz
    LINK=ftp://ftp.esrf.fr/tmp/spectrocrunch/$FILE

    # Download
    cd $CACHED_FOLDER
    if [ ! -f $FILE ]; then
        echo "Looking for pre-build libraries ..."
        NOSUCHFILE="$( wget --spider $LINK 2>&1 > /dev/null | grep No | wc -l)"
        if [ $NOSUCHFILE -eq 0 ]; then
            wget $LINK
        fi
    fi

    # Unpack
    if [ -f $FILE ]; then
        echo "Unpack pre-build libraries ..."
        tar -xzf $FILE -C $BUILD_FOLDER
    fi
    cd $BUILD_FOLDER
fi

