# 
# Install project dependencies (system and pypi).
# 

. $PSScriptRoot\funcs.ps1
. $PSScriptRoot\funcs-python.ps1
#. $PSScriptRoot\funcs-simpleelastix.ps1
. $PSScriptRoot\funcs-xraylib.ps1
#. $PSScriptRoot\funcs-pytmm.ps1
#. $PSScriptRoot\funcs-fdmnes.ps1


function install_system_dependencies()
{
    cprintstart
    cprint "Installing system requirements ..."

    if (!(dryrun)) {
        require_web_access
        pip_install numpy # silx/pyopencl

        # packages from Christoph Gohlke
        #python_init_compiler
        pip_install pipwin
        python_exec -m pipwin refresh
        python_exec -m pipwin install pyopencl # silx
        python_exec -m pipwin install shapely # spectrocrunch
    }
}


function install_system_dependencies_dev()
{
    cprintstart
    cprint "Installing system requirements (dev) ..."

    if (!(dryrun)) {
    }
}


function install_nopypi_dependencies()
{
    #require_simpleelastix
    #require_xraylib
    #require_pytmm
    #require_fdmnes
}


function install_pypi_dependencies()
{
    cprintstart
    cprint "Installing pypi requirements ..."
    if (!(dryrun)) {
        require_web_access
        pip_install "-r $(project_folder)\requirements.txt"
        pip_install pymca
    }
}


function install_pypi_dependencies_dev()
{
    cprintstart
    cprint "Installing pypi requirements (dev) ..."
    if (!(dryrun)) {
        require_web_access
        pip_install "-r $(project_folder)\requirements-dev.txt"
    }
}
