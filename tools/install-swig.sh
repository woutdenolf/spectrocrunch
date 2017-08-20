#!/bin/bash
# 
# Install swig on Linux.
# 

# ============Initialize environment============
SCRIPT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $SCRIPT_ROOT/funcs.sh
initEnv

# ============Install dependencies============
cprint "Install swig dependencies ..."
if [[ $NOTDRY == true && $SYSTEM_PRIVILIGES == true ]]; then
    mexec "apt-get -y install libpcre3 libpcre3-dev"
    mexec "apt-get -y remove swig2.0" # Otherwise you get conflicts
fi

# ============Install swig============
cprint "Download swig ..."
mkdir -p swig
cd swig
if [[ $NOTDRY == true && ! -d swig-3.0.12 ]]; then
    wget https://sourceforge.net/projects/swig/files/swig/swig-3.0.12/swig-3.0.12.tar.gz
    tar -xzf swig-3.0.12.tar.gz
fi
cd swig-3.0.12

if [[ ! -f Makefile ]]; then
    cprint "Configure swig ..."
    if [[ $NOTDRY == true ]]; then
        if [[ $INSTALL_SYSTEMWIDE == true ]]; then
            ./configure
        else
            mkdir -p $SPECTROCRUNCHLOCAL
            ./configure --prefix=$SPECTROCRUNCHLOCAL
        fi
    fi
    cprint "Build swig ..."
    if [[ $NOTDRY == true ]]; then
        make -s -j2
    fi
fi


cprint "Install swig ..."
ls
if [[ $NOTDRY == true ]]; then
    mmakeinstall

    # Add path just for this installation script
    addBinPath $SPECTROCRUNCHLOCAL/bin
fi

BUILDSTEP=$(( $BUILDSTEP+1 ))
BUILDSTEPS=$(( $BUILDSTEPS+1 ))
cd $INSTALL_WD

