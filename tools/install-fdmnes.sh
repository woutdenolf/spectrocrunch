#!/bin/bash
# 
# Install fdmnes on Linux.
# 

# ============Initialize environment============
if [ -z $PYTHONBIN ]; then
    PYTHONBIN=python
fi

if [ -z $NOTDRY ]; then
    NOTDRY=true
fi

if [ -z $BUILDSTEP ]; then
    BUILDSTEP=0
    BUILDSTEPS=0
fi

if [ -z $SYSTEM_PRIVILIGES ]; then
    if [[ -z "$((sudo -n true) 2>&1)" ]]; then
        export SYSTEM_PRIVILIGES=true 
    else
        export SYSTEM_PRIVILIGES=false
    fi
fi

# ============Install fdmnes============
mkdir -p fdmnes
cd fdmnes

FDMNESVERSION=`curl -s https://api.github.com/repos/gudasergey/fdmnes/releases | grep tag_name | head -n 1 | cut -d '"' -f 4`
FDMNESZIPNAME=fdmnes_${FDMNESVERSION}.zip
FDMNESLINK=http://neel.cnrs.fr/IMG/zip/${FDMNESZIPNAME}

echo -e "${hcol}Install fdmnes ...${ncol}"
if [[ $NOTDRY == true ]]; then
    if [[ $SYSTEM_PRIVILIGES == true ]]; then
        if [ -d /sware/exp/fdmnes ]; then
            # Use ESRF share
            if [ ! -d /opt/fdmnes ]; then
                ln -s /sware/exp/fdmnes /opt/fdmnes
            fi
            if [ ! -d /usr/local/bin/fdmnes ]; then
                ln -s /opt/fdmnes /usr/local/bin/fdmnes
            fi
        else
            # Install system-wide
            if [[ ! -f /usr/local/bin/fdmnes ]]; then
                curl -O ${FDMNESLINK}
                sudo -E unzip ${FDMNESZIPNAME} -d /usr/local
                sudo -E ln -s /usr/local/fdmnes/fdmnes_linux64 /usr/local/bin/fdmnes
            fi
        fi
    else
        # Install locally
        if [[ ! -f $HOME/.local/bin/fdmnes ]]; then
            curl -O ${FDMNESLINK}
            mkdir -p $HOME/.local/bin
            unzip ${FDMNESZIPNAME} -d $HOME/.local
            ln -s $HOME/.local/fdmnes/fdmnes_linux64 $HOME/.local/bin/fdmnes
        fi
    fi
fi

echo -e "${hcol}Download pyfdmnes ...${ncol}"
if [[ $NOTDRY == true && ! -d pyFDMNES ]]; then
    git clone https://github.com/woutdenolf/pyFDMNES.git pyFDMNES
fi

echo -e "${hcol}Install pyfdmnes ...${ncol}"
if [[ $NOTDRY == true ]]; then
    cd pyFDMNES
    $PYTHONBIN setup.py build -f
    $PYTHONBIN setup.py install
fi

BUILDSTEP=$(( $BUILDSTEP+1 ))
BUILDSTEPS=$(( $BUILDSTEPS+1 ))
