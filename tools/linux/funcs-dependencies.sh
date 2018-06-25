#!/bin/bash
# 
# Install project dependencies (system and pypi).
# 

SCRIPT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${SCRIPT_ROOT}/funcs.sh
source ${SCRIPT_ROOT}/funcs-python.sh
source ${SCRIPT_ROOT}/funcs-opencl.sh
source ${SCRIPT_ROOT}/funcs-opengl.sh
source ${SCRIPT_ROOT}/funcs-simpleelastix.sh
source ${SCRIPT_ROOT}/funcs-xraylib.sh
source ${SCRIPT_ROOT}/funcs-pytmm.sh
source ${SCRIPT_ROOT}/funcs-fdmnes.sh

function install_system_dependencies()
{
    cprintstart
    cprint "Installing system requirements ..."

    if [[ $(dryrun) == false ]]; then
        require_web_access

        pip_install numpy # silx

        pip_install cython # xraylib

        mapt-get libhdf5-serial-dev libhdf5-dev # h5py

        require_pyopencl # silx

        require_pyopengl # pymca

        require_pyqt # pymca

        #mapt-get libgl1-mesa-dev libglu1-mesa-dev mesa-common-dev # pymca
        #pip_install --egg pymca #TODO: wait for pymca to get fixed

        #mapt-get libgeos-dev # shapely
    fi
}


function install_system_dependencies_dev()
{
    cprintstart
    cprint "Installing system requirements (dev) ..."

    if [[ $(dryrun) == false ]]; then
        require_web_access

        mapt-get pandoc # nbsphinx
    fi
}


function install_pypi_dependencies()
{
    cprintstart
    cprint "Installing pypi requirements ..."
    if [[ $(dryrun) == false ]]; then
        require_web_access
        pip_install -r $(project_folder)/requirements.txt
    fi
}


function install_pypi_dependencies_dev()
{
    cprintstart
    cprint "Installing pypi requirements (dev) ..."
    if [[ $(dryrun) == false ]]; then
        require_web_access
        pip_install -r $(project_folder)/requirements-dev.txt
    fi
}


function install_nopypi_dependencies()
{
    require_simpleelastix
    require_xraylib
    require_pytmm
    require_fdmnes
}


