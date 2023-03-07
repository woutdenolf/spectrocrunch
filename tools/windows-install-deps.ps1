# 
# This script will install all Python and system dependencies.
# 

# ============Parse arguments============
param(
    [Parameter(Mandatory=$false)]
    [switch]$h = $false,
    [Parameter(Mandatory=$false)]
    [string]$v = "",
    [Parameter(Mandatory=$false)]
    [switch]$u = $false,
    [Parameter(Mandatory=$false)]
    [switch]$y = $false,
    [Parameter(Mandatory=$false)]
    [switch]$d = $false,
    [Parameter(Mandatory=$false)]
    [switch]$t = $false,
    [Parameter(Mandatory=$false)]
    [switch]$x = $false,
    [Parameter(Mandatory=$false)]
    [switch]$s = $false,
    [Parameter(Mandatory=$false)]
    [ValidateSet(0,32,64)]
    [int]$arch = 0
)

$global:ErrorActionPreference = "Stop"
$global:ARG_SYSTEMWIDE = !$u
$global:ARG_DRY = $d
$global:ARG_BASHRC = !$t
$global:ARG_FORCECHOICE = $y
$global:ARG_PYTHONV = $v
$global:ARG_DEV = $x
$global:ARG_ARCH = $arch
$global:ARG_SKIPLONG = $s
$global:ARG_HELP = $h

# ============Usage============
function show_help()
{
  echo "
        Usage: windows-install-deps  [-v version] [-y] [-d] [-x]

        -v version      Python version to be used (2, 3, 2.7, 3.5, ...).
        -y              Answer yes to everything.
        -d              Dry run.
        -u              Install for user only.
        -x              Dev install.
        -t              Do not modify bashrc.
        -s              Skip libraries that take long to compile.

        For Example: ./windows-install-deps.sh -v 3 -d

        -h              Help
        
       Environment variables(optional):
       
        PROJECT_PREFIX: directory to install dependencies
        PROJECT_RESOURCE_DIR: directory of the project resource file
       "
}


function main()
{
    if ($global:ARG_HELP) {
        show_help
        return
    }
    
    # ============Initialize environment============
    $local:GLOBAL_WD = Get-Location
    . $PSScriptRoot\windows\funcs-python.ps1
    dryrun reset $global:ARG_DRY
    install_systemwide reset $global:ARG_SYSTEMWIDE
    install_arch reset $global:ARG_ARCH
    if ((dryrun)) {
        cprint "This is a dry run."
    }

    # ============Ask for confirmation to proceed============
    if (!(dryrun)) {
        require_python $global:ARG_PYTHONV
        require_pip
    }
    install_arch reset $(python_arch)

    python_info 
    pip_info
    install_info

    $local:choice = $true
    if (!$global:ARG_FORCECHOICE) {
        $local:choice = (YesNoQuestion "Continue installation?")
    }
    if (!$local:choice) {
        return
    }
    timer reset

    # ============Install dependencies============
    . $PSScriptRoot\windows\funcs-dependencies

    mkdir $(python_depdir) -Force
    cd $(python_depdir)

    install_system_dependencies
    install_pypi_dependencies
    install_nopypi_dependencies

    if ($global:ARG_DEV) {
        install_system_dependencies_dev
        install_pypi_dependencies_dev
    }

    # ============Cleanup============
    cprintstart
    cprint "Cleaning up ..."

    if ((dryrun)) {
        cprint "This was a dryrun, so it didn't install anything."
    } else {
        cprint "All done! You should now be able to install $(project_name)."
    }

    timer
    cd $local:GLOBAL_WD
}


main