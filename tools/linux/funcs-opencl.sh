#!/bin/bash
# 
# Install pyopencl and drivers.
# 

SCRIPT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${SCRIPT_ROOT}/funcs.sh
source ${SCRIPT_ROOT}/funcs-python.sh


function intel_name()
{
    echo "opencl_runtime_16.1.2_x64_rh_6.4.0.37"
}

function amd_name()
{
    echo "AMD-APP-SDKInstaller-v3.0.130.136-GA-linux64"
}


function intel_url()
{
    # https://software.intel.com/en-us/articles/opencl-drivers#latest_CPU_runtime
    echo "http://registrationcenter-download.intel.com/akdlm/irc_nas/12556/$(intel_name).tgz"
}

function amd_url()
{
    echo "http://debian.nullivex.com/amd/$(amd_name).tar.bz2"
}


function intel_build_dependencies()
{
    mapt-get install cpio
}


function intel_install_opencldrivers()
{
    local restorewd=$(pwd)

    cprint "Download Intel OpenCL runtime ..."
    mkdir -p opencl
    cd opencl

    local PACKAGE_NAME=$(intel_name)
    if [[ $(dryrun) == false && ! -d ${PACKAGE_NAME} ]]; then
        require_web_access
        curl -L $(intel_url) --output ${PACKAGE_NAME}.tar.gz
        mkdir -p ${PACKAGE_NAME}
        tar -xzf ${PACKAGE_NAME}.tar.gz -C ${PACKAGE_NAME} --strip-components=1
        rm -f ${PACKAGE_NAME}.tar.gz
    fi

    local prefix=$(project_opt)/opencl
    local prefixstr=$(project_optstr)/opencl

    cprint "Install Intel OpenCL runtime ..."
    if [[ $(dryrun) == false ]]; then
        cd ${PACKAGE_NAME}
        intel_build_dependencies

        sed 's/decline/accept/g' -i silent.cfg
        local repl=$(strreplace ${prefix} "/" "\/")
        sed "s/=\/opt/=${repl}/g" -i silent.cfg
        sed 's/DEFAULTS/ALL/g' -i silent.cfg

        # in subshell to capture the exit statements
        if [[ $(system_privileges) == true ]]; then
            (sudo -E ./install.sh -s silent.cfg) 
        else
            (./install.sh -s silent.cfg) # this will not work
        fi
        
        addProfile $(project_resource) "# Installed Intel OpenCL runtime: ${prefixstr}"
        addLibPath ${prefix}/lib
        addLibPathProfile $(project_resource) "${prefixstr}/intel/opencl/lib64"
    fi

    cd ${restorewd}
}


function amd_install_opencldrivers()
{
    local restorewd=$(pwd)

    cprint "Download AMD OpenCL SDK ..."
    mkdir -p opencl
    cd opencl

    local PACKAGE_NAME=$(amd_name)
    if [[ $(dryrun) == false && ! -d ${PACKAGE_NAME} ]]; then
        require_web_access
        curl -L $(amd_url) --output ${PACKAGE_NAME}.tar.bz2
        mkdir -p ${PACKAGE_NAME}
        tar -xjf ${PACKAGE_NAME}.tar.bz2 -C ${PACKAGE_NAME}
        rm -f ${PACKAGE_NAME}.tar.bz2
    fi

    local prefix=$(project_opt)/opencl/AMDAPPSDK
    local prefixstr=$(project_optstr)/opencl/AMDAPPSDK

    cprint "Install AMD OpenCL SDK ..."
    if [[ $(dryrun) == false ]]; then
        cd ${PACKAGE_NAME}
        
        # Install in opt
        mkdir -p ${prefix}
        sh AMD-APP-SDK*.sh --tar -xf -C ${prefix}

        # Register the ICD
        local OPENCL_VENDOR_PATHSTR
        export OPENCL_VENDOR_PATH=${prefix}/etc/OpenCL/vendors
        OPENCL_VENDOR_PATHSTR="${prefixstr}/etc/OpenCL/vendors"
        mkdir -p ${OPENCL_VENDOR_PATH}
        if [[ $(os_arch) == 64 ]]; then
            if [ ! -f ${OPENCL_VENDOR_PATH}/amdocl64.icd ]; then
                echo libamdocl64.so > ${OPENCL_VENDOR_PATH}/amdocl64.icd
            fi
        else
            if [ ! -f ${OPENCL_VENDOR_PATH}/amdocl32.icd ]; then
                echo libamdocl32.so > ${OPENCL_VENDOR_PATH}/amdocl32.icd
            fi
        fi

        # Environment variables
        addProfile $(project_resource) "# Installed AMD SDK: ${prefixstr}"
        addProfile $(project_resource) "export OPENCL_VENDOR_PATH=${OPENCL_VENDOR_PATHSTR}"
        if [[ $(os_arch) == 64 ]]; then
            addLibPath ${prefix}/lib/x86_64
            addLibPathProfile $(project_resource) "${prefixstr}/lib/x86_64"
        else
            addLibPath ${prefix}/lib/x86
            addLibPathProfile $(project_resource) "${prefixstr}/lib/x86"
        fi
        addInclPath ${prefix}/include

        # Show info
        if [[ $(os_arch) == 64 ]]; then
            chmod +x ${prefix}/bin/x86_64/clinfo
            ${prefix}/bin/x86_64/clinfo
        else
            chmod +x ${prefix}/bin/x86/clinfo
            ${prefix}/bin/x86/clinfo
        fi

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
        mapt-get install $(python_apt)-mako
        mapt-get install $(python_apt)-pyopencl
        if [[ $(pyopencl_test) == true ]]; then
            return
        fi

        # ICD = “installable client driver”
        # ocl-icd-libopencl1: ICD loader
        # ocl-icd-opencl-dev: ICD loader dev
        # opencl-headers: C/C++ headers for OpenCL API
        # libffi-dev: 
        # mako: templating language for python
        # pybind11: C++/Python connection
        
        mapt-get install ocl-icd-opencl-dev ocl-icd-libopencl1 opencl-headers libffi-dev
        
        intel_install_opencldrivers
        pip_install numpy pybind11 mako
        pip_install pyopencl

        if [[ $(pyopencl_test) == false ]]; then
            amd_install_opencldrivers
            pip_install numpy pybind11 mako
            pip_install pyopencl
        fi
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
