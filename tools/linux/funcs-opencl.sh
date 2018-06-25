#!/bin/bash
# 
# Install pyopencl and drivers.
# 

SCRIPT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${SCRIPT_ROOT}/funcs.sh
source ${SCRIPT_ROOT}/funcs-python.sh

function intel_url()
{
    # https://software.intel.com/en-us/articles/opencl-drivers#latest_CPU_runtime
    echo "http://registrationcenter-download.intel.com/akdlm/irc_nas/12556/opencl_runtime_16.1.2_x64_rh_6.4.0.37.tgz"
}


function intel_build_dependencies()
{
    mapt-get cpio
}


function intel_install_opencldrivers()
{
    local restorewd=$(pwd)

    cprint "Download OpenCL runtime ..."
    mkdir -p opencl
    cd opencl

    local PACKAGE_NAME=opencl_runtime_16.1.2_x64_rh_6.4.0.37
    if [[ $(dryrun) == false  && ! -d${PACKAGE_NAME} ]]; then
        require_web_access
        wget -q $(intel_url) -O opencl_runtime.tgz
        tar -xzf opencl_runtime.tgz
    fi

    cd ${PACKAGE_NAME}

    local prefix=$(project_opt)/opencl
    local prefixstr=$(project_optstr)/opencl

    cprint "Install OpenCL runtime ..."
    if [[ $(dryrun) == false ]]; then
        intel_build_dependencies

        sed 's/decline/accept/g' -i silent.cfg
        local repl=$(strreplace ${prefix} "/" "\/")
        sed "s/=\/opt/=${repl}/g" -i silent.cfg
        sed 's/DEFAULTS/ALL/g' -i silent.cfg

        (sudo -E ./install.sh -s silent.cfg) # in subshell to capture the exit statements

        addProfile $(project_resource) "# Installed xraylib: ${prefixstr}"
        addLibPath ${prefix}/lib
        addLibPathProfile $(project_resource) "${prefixstr}/intel/opencl/lib64"
    fi

    cd ${restorewd}
}


function pyopencl_test()
{
    python_get $'try:\n import pyopencl\n assert pyopencl.get_platforms()\n print("true")\nexcept Exception:\n print("false")'
}


function pyopencl_install()
{
    if [[ $(dryrun) == false ]]; then
        mapt-get $(python_bin)-pyopencl
        if [[ $(pyopencl_test) == true ]]; then
            return
        fi

        pip_install pyopencl
        if [[ $(pyopencl_test) == true ]]; then
            return
        fi

        mapt-get ocl-icd-libopencl1 opencl-headers ocl-icd-opencl-dev libffi-dev
        intel_install_opencldrivers
    fi
}


function require_pyopencl()
{
    cprintstart
    cprint "Verify pyopencl ..."

    # Requirements (for running)
    require_python

    # Check
    if [[ $(pyopencl_test) == true ]]; then
        cprint "Python module \"pyopencl\" is working"
        cprintend
        return
    fi

    # Install
    pyopencl_install

    # Check
    if [[ $(pyopencl_test) == true ]]; then
        cprint "Python module \"pyopencl\" is working"
    else
        cprint "Python module \"pyopencl\" is NOT working"
    fi

    cprintend
}
