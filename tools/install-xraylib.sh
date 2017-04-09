#!/bin/bash
# 
# Install xraylib on Linux.
# 

# ============Initialize environment============
SCRIPT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $SCRIPT_ROOT/funcs.sh
initEnv

# ============Install dependencies============
cprint "Install xraylib dependencies ..."
if [[ $NOTDRY == true && $SYSTEM_PRIVILIGES == true ]]; then
    sudo -E apt-get -y install swig
    sudo -E apt-get -y install $PYTHONBINAPT-dev
fi

# ============Install xraylib============
if [ ! -f xraylib/xraylib-3.2.0/python/.libs/_xraylib.so ]; then
    cprint "Download xraylib ..."
    mkdir -p xraylib
    cd xraylib
    if [[ $NOTDRY == true && ! -d xraylib-3.2.0 ]]; then
        curl -O http://lvserver.ugent.be/xraylib/xraylib-3.2.0.tar.gz
        tar -xvf xraylib-3.2.0.tar.gz
        cd xraylib-3.2.0
    fi

    cprint "Configure xraylib ..."
    if [[ $NOTDRY == true ]]; then
        if [[ $SYSTEM_PRIVILIGES == true ]]; then
            ./configure --enable-python \
                        --enable-python-integration \
                        --disable-java \
                        --disable-lua \
                        --disable-ruby \
                        --disable-php \
                        --disable-pascal \
                        --disable-idl \
                        --disable-perl \
                        PYTHON=$PYTHON_EXECUTABLE \
                        PYTHON_VERSION=$PYTHONV
        else
            mkdir -p $SPECTROCRUNCHLOCAL
            ./configure --prefix="$SPECTROCRUNCHLOCAL" \
                        --enable-python \
                        --enable-python-integration \
                        --disable-java \
                        --disable-lua \
                        --disable-ruby \
                        --disable-php \
                        --disable-pascal \
                        --disable-idl \
                        --disable-perl \
                        PYTHON=$PYTHON_EXECUTABLE \
                        PYTHON_VERSION=$PYTHONV
        fi
    fi

    cprint "Build xraylib ..."
    if [[ $NOTDRY == true ]]; then
        make -s -j2
        make check
    fi
else
    cd xraylib/xraylib-3.2.0
fi

cprint "Install xraylib ..."
if [[ $NOTDRY == true ]]; then
    mexec "make install -s"

    addProfile $SPECTROCRUNCHRC "# Installed xraylib: $SPECTROCRUNCHOPTSTR/xraylib"
    addLibPath $SPECTROCRUNCHLOCAL/lib
    addLibPathProfile $SPECTROCRUNCHRC "$SPECTROCRUNCHLOCALSTR/lib"
    addBinPath $SPECTROCRUNCHLOCAL/bin
    addBinPathProfile $SPECTROCRUNCHRC "$SPECTROCRUNCHLOCALSTR/bin"
fi

BUILDSTEP=$(( $BUILDSTEP+1 ))
BUILDSTEPS=$(( $BUILDSTEPS+1 ))
cd $INSTALL_WD

