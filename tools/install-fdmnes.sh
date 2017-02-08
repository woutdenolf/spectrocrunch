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
                curl -O http://neel.cnrs.fr/IMG/zip/fdmnes_2017_01_10.zip
                sudo -E unzip fdmnes_2017_01_10.zip -d /usr/local
                sudo -E ln -s /usr/local/fdmnes/fdmnes_linux64 /usr/local/bin/fdmnes
            fi
        fi
    else
        # Install locally
        if [[ ! -f $HOME/.local/bin/fdmnes ]]; then
            curl -O http://neel.cnrs.fr/IMG/zip/fdmnes_2017_01_10.zip
            mkdir -p $HOME/.local/bin
            unzip fdmnes_2017_01_10.zip -d $HOME/.local
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

