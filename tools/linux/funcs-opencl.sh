#!/bin/bash
# 
# Install pyopencl and drivers.
# 

SCRIPT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${SCRIPT_ROOT}/funcs.sh
source ${SCRIPT_ROOT}/funcs-python.sh


function amd_name()
{
    echo "AMD-APP-SDKInstaller-v3.0.130.136-GA-linux64"
}


function intel_url()
{
    # https://software.intel.com/en-us/articles/opencl-drivers#latest_CPU_runtime
    echo "https://registrationcenter-download.intel.com/akdlm/irc_nas/vcp/15532/$(intel_name).tgz"
}


function amd_url()
{
    echo "http://cs.wisc.edu/~riccardo/assets/$(amd_name).tar.bz2"
}


function intel_install_opencldrivers()
{
    local restorewd=$(pwd)

    cprint "Download Intel OpenCL runtime ..."
    mkdir -p opencl
    cd opencl

    local PACKAGE_NAME="intel-opencl_21.01.18793_amd64"
    if [[ $(dryrun) == false && ! -d ${PACKAGE_NAME} ]]; then
        require_web_access
        wget https://github.com/intel/compute-runtime/releases/download/21.01.18793/${PACKAGE_NAME}.deb
        dpkg-deb -xv ${PACKAGE_NAME}.deb ${PACKAGE_NAME}
    fi

    local PACKAGE_NAME2="intel-gmmlib_20.3.2_amd64"
    if [[ $(dryrun) == false && ! -d ${PACKAGE_NAME2} ]]; then
        require_web_access
        wget https://github.com/intel/compute-runtime/releases/download/21.01.18793/${PACKAGE_NAME2}.deb
        dpkg-deb -xv ${PACKAGE_NAME2}.deb ${PACKAGE_NAME2}
    fi

    local prefix=$(project_opt)/opencl
    local prefixstr=$(project_optstr)/opencl

    cprint "Install Intel OpenCL runtime ..."
    if [[ $(dryrun) == false ]]; then
        # Register the ICD
        export OCL_ICD_VENDORS=${prefix}/etc/OpenCL/vendors
        local OCL_ICD_VENDORS_STR="${prefixstr}/etc/OpenCL/vendors"
        mkdir -p ${OCL_ICD_VENDORS}
        echo libigdrcl.so > ${OCL_ICD_VENDORS}/intel.icd

        # Install the libraries
        mkdir -p ${prefix}/lib/x86_64
        cp -a ${PACKAGE_NAME}/usr/local/lib/intel-opencl/* ${prefix}/lib/x86_64
        cp -a ${PACKAGE_NAME2}/usr/local/lib/* ${prefix}/lib/x86_64
        chmod 755 ${prefix}/lib/x86_64/*

        # Environment variables
        addProfile $(project_resource) "# Installed Intel OpenCL runtime: ${prefixstr}"
        addProfile $(project_resource) "export OCL_ICD_VENDORS=${OCL_ICD_VENDORS_STR}"
        addLibPath "${prefix}/lib/x86_64"
        addLibPathProfile $(project_resource) "${prefixstr}/lib/x86_64"
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
        export OCL_ICD_VENDORS=${prefix}/etc/OpenCL/vendors
        local OCL_ICD_VENDORS_STR="${prefixstr}/etc/OpenCL/vendors"
        mkdir -p ${OCL_ICD_VENDORS}
        if [[ $(os_arch) == 64 ]]; then
            if [ ! -f ${OCL_ICD_VENDORS}/amdocl64.icd ]; then
                echo libamdocl64.so > ${OCL_ICD_VENDORS}/amdocl64.icd
            fi
        else
            if [ ! -f ${OCL_ICD_VENDORS}/amdocl32.icd ]; then
                echo libamdocl32.so > ${OCL_ICD_VENDORS}/amdocl32.icd
            fi
        fi

        # Environment variables
        addProfile $(project_resource) "# Installed AMD SDK: ${prefixstr}"
        addProfile $(project_resource) "export OCL_ICD_VENDORS=${OCL_ICD_VENDORS_STR}"
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

        if [[ -d /etc/OpenCL/vendors/ ]];then
            # Needs to be done explicitely on slurm, not sure why
            export OCL_ICD_VENDORS=/etc/OpenCL/vendors/
            if [[ $(pyopencl_test) == true ]]; then
                addProfile $(project_resource) "# Use system installed OpenCL ICD's"
                addProfile $(project_resource) "export OCL_ICD_VENDORS=${OCL_ICD_VENDORS}"
                return
            fi
            unset OCL_ICD_VENDORS
        fi

        # ICD = “installable client driver”
        # ocl-icd-libopencl1: ICD loader
        # ocl-icd-opencl-dev: ICD loader dev
        # opencl-headers: C/C++ headers for OpenCL API
        # libffi-dev: 
        # mako: templating language for python
        # pybind11: C++/Python connection
        
        mapt-get install ocl-icd-opencl-dev
        mapt-get install ocl-icd-libopencl1
        mapt-get install opencl-headers
        mapt-get install libffi-dev

        intel_install_opencldrivers
        pip_install numpy
        pip_install pybind11
        pip_install mako
        pip_install pyopencl

        return # AMD no longer supported

        if [[ $(pyopencl_test) == false ]]; then
            amd_install_opencldrivers
            pip_install numpy
            pip_install pybind11
            pip_install mako
            pip_install pyopencl
        fi
    fi
}


function require_pyopencl()
{
    cprintstart "Require pyopencl"

    # Requirements (for running)
    require_python

    # Check
    if [[ $(pyopencl_test) == true ]]; then
        cprint "Python module \"pyopencl\" is working"
        cprintend "Require pyopencl"
        return
    fi

    # Install
    pyopencl_install

    # Check
    if [[ $(pyopencl_test) == true ]]; then
        cprint "Python module \"pyopencl\" is working"
    else
        cerror "Python module \"pyopencl\" is NOT working"
    fi

    cprintend "Require pyopencl"
}
