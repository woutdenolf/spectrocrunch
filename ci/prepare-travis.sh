#!/bin/bash
# 
# This script will install all spectrocrunch Python 2 and 3 dependencies for CI.
# 

# Download pre-build libraries
PYTHONV=`python -c "import sys;t='{v[0]}.{v[1]}'.format(v=list(sys.version_info[:2]));print(t)";`

if [ ! -d ${PYTHONV} ]; then
    FILE=spectrocrunch.travis.python${PYTHONV}.tgz
    if [[ ${PYTHONV} == "2.7" ]]; then
        LINK=https://transfer.sh/fGRlO/$FILE
    else
        http://ftp.esrf.fr/tmp/$FILE
    fi

    # Download
    cd $CACHED_FOLDER
    if [ ! -f $FILE ]; then
        echo "Download pre-build libraries ..."
        wget ${LINK}
    fi

    # Unpack
    if [ -f $FILE ]; then
        echo "Unpack pre-build libraries ..."
        tar -xzf $FILE -C $BUILD_FOLDER
    fi
    cd $BUILD_FOLDER
fi

