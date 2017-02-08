#!/bin/bash
# 
# Install xraylib on Linux.
# 

# ============Initialize environment============
if [ -z $NOTDRY ]; then
    NOTDRY=true
fi

if [ -z $BUILDSTEP ]; then
    BUILDSTEP=0
fi

if [ -z $TIMELEFT ]; then
    TIMELEFT=true
fi

if [ -z $TIMELIMITED ]; then
    TIMELIMITED=false
fi

if [ -z $SYSTEM_PRIVILIGES ]; then
    if [[ -z "$((sudo -n true) 2>&1)" ]]; then
        export SYSTEM_PRIVILIGES=true 
    else
        export SYSTEM_PRIVILIGES=false
    fi
fi

# ============Install cmake============
echo -e "${hcol}Download cmake ...${ncol}"
mkdir -p cmake
cd cmake
if [[ $NOTDRY == true && ! -d cmake-3.7.2 ]]; then
    wget --no-check-certificate http://www.cmake.org/files/v3.7/cmake-3.7.2.tar.gz
    tar -xzf cmake-3.7.2.tar.gz
fi
cd cmake-3.7.2

if [[ $TIMELEFT == true && ! -f Makefile ]]; then
    echo -e "${hcol}Configure cmake ...${ncol}"
    if [[ $NOTDRY == true ]]; then
        if [[ $SYSTEM_PRIVILIGES == true ]]; then
            ./configure
        else
            mkdir -p $HOME/.local
            ./configure --prefix=$HOME/.local
        fi
    fi
    echo -e "${hcol}Build cmake ...${ncol}"
    if [[ $NOTDRY == true ]]; then
        make -s -j2
    fi

    if [[ $TIMELIMITED == true ]]; then
        TIMELEFT=false
    fi
fi

echo -e "${hcol}Install cmake ...${ncol}"
if [[ $NOTDRY == true ]]; then
    if [[ $SYSTEM_PRIVILIGES == true ]]; then
        sudo -E make install -s
    else
        make install -s
    fi
fi

BUILDSTEP=$(( $BUILDSTEP+1 ))

