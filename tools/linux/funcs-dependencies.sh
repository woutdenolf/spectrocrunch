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
source ${SCRIPT_ROOT}/funcs-xmimsim.sh
source ${SCRIPT_ROOT}/funcs-xrmc.sh
source ${SCRIPT_ROOT}/funcs-pytmm.sh
source ${SCRIPT_ROOT}/funcs-fdmnes.sh

function install_system_dependencies()
{
    cprintstart "Installing system requirements"
    if [[ $(dryrun) == false ]]; then
        require_web_access
        pip_install numpy # silx
        pip_install cython # xraylib
        mapt-get install libhdf5-serial-dev libhdf5-dev # h5py
        mapt-get install $(python_apt)-tk # matplotlib
        require_pyopencl # silx
        require_pyopengl # pymca
        require_pyqt # pymca
        #mapt-get install libgeos-dev # shapely
    fi
    cprintend "Installing system requirements"
}


function install_system_dependencies_dev()
{
    cprintstart "Installing system requirements (dev)"
    if [[ $(dryrun) == false ]]; then
        require_web_access
        mapt-get install pandoc # nbsphinx
    fi
    cprint "Finish: Installing system requirements (dev)"
    cprintend "Installing system requirements (dev)"
}


function install_nopypi_dependencies()
{
    require_xraylib
    require_xmimsim
    require_xrmc
    require_simpleelastix
    require_pytmm
    if [[ $(python_major_version) == 2 ]];then
        require_fdmnes
    fi
}


function install_pypi_dependencies()
{
    cprintstart "Installing pypi requirements"
    if [[ $(dryrun) == false ]]; then
        require_web_access
        pip_install -r "$(project_folder)/requirements.txt"
        pip_install pymca
    fi
    cprintend "Installing pypi requirements"
}


function install_pypi_dependencies_dev()
{
    cprintstart "Installing pypi requirements (dev)"
    if [[ $(dryrun) == false ]]; then
        require_web_access
        pip_install -r "$(project_folder)/requirements-dev.txt"
    fi
    cprintend "Installing pypi requirements (dev)"
}

