#!/bin/bash
# 
# Install xraylib on Linux.
# 

# ============Initialize environment============
if [ -z $PYTHONBIN ]; then
    PYTHONBIN=python
fi
if [ -z $PYTHONBINAPT ]; then
    PYTHONBINAPT=python
fi
PYTHON_EXECUTABLE=$(which $PYTHONBIN) # full path
PYTHONV=`$PYTHONBIN -c "import sys;t='{v[0]}.{v[1]}'.format(v=list(sys.version_info[:2]));print(t)";` # e.g. 2.7

if [ -z $NOTDRY ]; then
    NOTDRY=true
fi

if [ -z $BUILDSTEP ]; then
    BUILDSTEP=0
fi

if [ -z $SYSTEM_PRIVILIGES ]; then
    if [[ -z "$((sudo -n true) 2>&1)" ]]; then
        export SYSTEM_PRIVILIGES=true 
    else
        export SYSTEM_PRIVILIGES=false
    fi
fi

# ============Install dependencies============
echo -e "${hcol}Install xraylib dependencies ...${ncol}"
if [[ $NOTDRY == true && $SYSTEM_PRIVILIGES == true ]]; then
    sudo -E apt-get -y install swig
    sudo -E apt-get install $PYTHONBINAPT-dev
fi

# ============Install xraylib============
if [ ! -f xraylib/xraylib-3.2.0/python/.libs/_xraylib.so ]; then
    echo -e "${hcol}Download xraylib ...${ncol}"
    mkdir -p xraylib
    cd xraylib
    if [[ $NOTDRY == true && ! -d xraylib-3.2.0 ]]; then
        curl -O http://lvserver.ugent.be/xraylib/xraylib-3.2.0.tar.gz
        tar -xvf xraylib-3.2.0.tar.gz
        cd xraylib-3.2.0
    fi

    echo -e "${hcol}Configure xraylib ...${ncol}"
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
            mkdir -p $HOME/.local
            ./configure --prefix="$HOME/.local" \
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

    echo -e "${hcol}Build xraylib ...${ncol}"
    if [[ $NOTDRY == true ]]; then
        make -s -j2
        make check
    fi
else
    cd xraylib/xraylib-3.2.0
fi

echo -e "${hcol}Install xraylib ...${ncol}"
if [[ $NOTDRY == true ]]; then
    if [[ $SYSTEM_PRIVILIGES == true ]]; then
        sudo -E make install -s
    else
        make install -s
    fi
fi

BUILDSTEP=$(( $BUILDSTEP+1 ))

