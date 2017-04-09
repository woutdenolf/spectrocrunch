#!/bin/bash
# 
# Install opencl on Linux.
# 

# ============Initialize environment============
SCRIPT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $SCRIPT_ROOT/funcs.sh
initEnv

# ============Install AMD SDK============
#export AMDAPPSDKROOT=$SPECTROCRUNCHOPT/AMDAPPSDK
#AMDAPPSDKROOTSTR=$SPECTROCRUNCHOPTSTR/AMDAPPSDK
#if [[ $SYSTEM_PRIVILIGES == true ]]; then    
#    export OPENCL_VENDOR_PATH=/etc/OpenCL/vendors
#    OPENCL_VENDOR_PATHSTR="/etc/OpenCL/vendors"
#else
#    export OPENCL_VENDOR_PATH=$AMDAPPSDKROOT/etc/OpenCL/vendors
#    OPENCL_VENDOR_PATHSTR="\$AMDAPPSDKROOT/etc/OpenCL/vendors"
#fi

# TODO: for Travis
export AMDAPPSDKROOT=$HOME/.local/AMDAPPSDK
AMDAPPSDKROOTSTR="\$HOME/.local/AMDAPPSDK"
export OPENCL_VENDOR_PATH=$AMDAPPSDKROOT/etc/OpenCL/vendors
OPENCL_VENDOR_PATHSTR="\$AMDAPPSDKROOT/etc/OpenCL/vendors"

mkdir -p amd
cd amd

if [ ! -f AMD-APP-SDK-v2.9-1.599.381-GA-linux64.sh ]; then
    cprint "Download AMD SDK for opencl ..."
    
    if [[ $NOTDRY == true ]]; then
        bash $SCRIPT_ROOT/amd_sdk.sh

        # Untar
        tar -xjf AMD-SDK.tar.bz2
	fi
fi

cprint "Install AMD SDK ..."
if [[ $NOTDRY == true ]]; then
    # Install in opt
    mkdir -p $AMDAPPSDKROOT
    sh AMD-APP-SDK*.sh --tar -xf -C $AMDAPPSDKROOT

    # Register the ICD
    mkdir -p $OPENCL_VENDOR_PATH
    if [ "$OS_ARCH" == "64" ]; then
        if [ ! -f $OPENCL_VENDOR_PATH/amdocl64.icd ]; then
            echo libamdocl64.so > $OPENCL_VENDOR_PATH/amdocl64.icd
        fi
    else
        if [ ! -f $OPENCL_VENDOR_PATH/amdocl32.icd ]; then
            echo libamdocl32.so > $OPENCL_VENDOR_PATH/amdocl32.icd
        fi
    fi

    # Environment variables
    addProfile $SPECTROCRUNCHRC "# Installed AMD SDK: $AMDAPPSDKROOTSTR"
    addProfile $SPECTROCRUNCHRC "export AMDAPPSDKROOT=$AMDAPPSDKROOTSTR"
    addProfile $SPECTROCRUNCHRC "export OPENCL_VENDOR_PATH=$OPENCL_VENDOR_PATHSTR"
    if [ "$OS_ARCH" == "64" ]; then
        addLibPath $AMDAPPSDKROOT/lib/x86_64
        addLibPathProfile $SPECTROCRUNCHRC "\$AMDAPPSDKROOT/lib/x86_64"
    else
        addLibPath $AMDAPPSDKROOT/lib/x86
        addLibPathProfile $SPECTROCRUNCHRC "\$AMDAPPSDKROOT/lib/x86"
    fi
    addInclPath $AMDAPPSDKROOT/include

    # Show info
    if [ "$OS_ARCH" == "64" ]; then
        chmod +x $AMDAPPSDKROOT/bin/x86_64/clinfo
        $AMDAPPSDKROOT/bin/x86_64/clinfo
    else
        chmod +x $AMDAPPSDKROOT/bin/x86/clinfo
        $AMDAPPSDKROOT/bin/x86/clinfo
    fi
fi

BUILDSTEP=$(( $BUILDSTEP+1 ))
BUILDSTEPS=$(( $BUILDSTEPS+1 ))
cd $INSTALL_WD

