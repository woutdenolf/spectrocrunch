#!/bin/bash
# 
# This script will install all spectrocrunch Python 2 and 3 dependencies for CI.
# 

# ============Initialize environment============
hcol='\033[0;35m'
ncol='\033[0m'

RESTORE_WD=$(pwd)
SPECTROCRUNCH_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"/..

# ============AMD SDK============
mkdir -p amd
cd amd

export AMDAPPSDK=$(pwd)/AMDAPPSDK
export OPENCL_VENDOR_PATH=$AMDAPPSDK/etc/OpenCL/vendors
export LD_LIBRARY_PATH=$AMDAPPSDK/lib/x86_64:$LD_LIBRARY_PATH

if [ ! -f $OPENCL_VENDOR_PATH/amdocl64.icd ]; then
  echo -e "${hcol}Download AMD SDK for opencl ...${ncol}"
  bash $SPECTROCRUNCH_ROOT/ci/amd_sdk.sh

  echo -e "${hcol}Install AMD SDK ...${ncol}"
  tar -xjf AMD-SDK.tar.bz2
  
  mkdir -p $AMDAPPSDK
  sh AMD-APP-SDK*.sh --tar -xf -C $AMDAPPSDK

  mkdir -p $OPENCL_VENDOR_PATH
  echo libamdocl64.so > $OPENCL_VENDOR_PATH/amdocl64.icd

  chmod +x $AMDAPPSDK/bin/x86_64/clinfo
  $AMDAPPSDK/bin/x86_64/clinfo
fi

# ============Cleanup============
#echo -e "${hcol}Cleaning up ...${ncol}"
cd $RESTORE_WD


