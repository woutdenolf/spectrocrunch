# 
# Install project dependencies (system and pypi).
# 

. $PSScriptRoot\funcs.ps1
. $PSScriptRoot\funcs-python.ps1


function install_system_dependencies()
{
    cprintstart
    cprint "Installing system requirements ..."

    if (!(dryrun)) {
        require_web_access
        python_init_compiler
        pip_install numpy # silx/pyopencl

        # packages from Christoph Gohlke
        # TODO: check for compiler
        pip_install pipwin
        pipwin install pyopencl # silx
        pipwin install shapely # spectrocrunch
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
    # require...
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
