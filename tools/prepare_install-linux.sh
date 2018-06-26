#!/bin/bash
# 
# This script will install all system and python dependencies.
# 

# ============Usage============
show_help()
{
  echo "
        Usage: prepare_installation  [-v version] [-y] [-d] [-x]

        -v version      Python version to be used (2, 3, 2.7, 3.5, ...).
        -y              Answer yes to everything.
        -d              Dry run.
        -u              Install for user only.
        -x              Dev install.

        For Example: ./prepare_installation -v 3 -d

        -h              Help
       "
}

# ============Parse script arguments============
GLOBAL_SCRIPT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
GLOBAL_WD=$(pwd)
OPTIND=0
ARG_SYSTEMWIDE=true
ARG_DRY=false
ARG_FORCECHOICE=false
ARG_PYTHONV=""
ARG_DEV=false
while getopts "v:uyhdx" opt; do
  case $opt in
    h)
      show_help
      return 0
      ;;
    y)
      ARG_FORCECHOICE=true
      ;;
    u)
      ARG_SYSTEMWIDE=false
      ;;
    d)
      ARG_DRY=true
      ;;
    x)
      ARG_DEV=true
      ;;
    v)
      ARG_PYTHONV=${OPTARG}
      ;;
    \?)
      echo "Invalid option: -${OPTARG}. Use -h flag for help." >&2
      return 1
      ;;
  esac
done

# ============Initialize environment============
source $GLOBAL_SCRIPT_ROOT/linux/funcs-python.sh
dryrun reset ${ARG_DRY}
if [[ $(dryrun) == true ]]; then
    cprint "This is a dry run."
fi
install_systemwide reset ${ARG_SYSTEMWIDE}

# ============Ask for confirmation to proceed============
if [[ $(dryrun) == false ]]; then
    require_python ${ARG_PYTHONV}
    require_pip
fi

python_info 
pip_info
cprint "Root priviliges: $(system_privileges)"
cprint "System wide installation: $(install_systemwide)"
cprint "Prefix for dependencies: $(project_prefix)"
cprint "Opt directory: $(project_opt)"

if [[ ${ARG_FORCECHOICE} == false ]]; then
    read -p "Continue (Y/n)?" ARG_CHOICE
else
    ARG_CHOICE=${ARG_FORCECHOICE}
fi
case "${ARG_CHOICE}" in 
  y|Y ) ;;
  n|N ) 
        return 1;;
  * ) ;;
esac
timer reset

# ============Install dependencies============
source $GLOBAL_SCRIPT_ROOT/linux/funcs-dependencies.sh

mkdir -p dep_$(python_full_version)
cd dep_$(python_full_version)

install_system_dependencies
install_pypi_dependencies
install_nopypi_dependencies

if [[ ${ARG_DEV} == true ]]; then
    install_system_dependencies_dev
    install_pypi_dependencies_dev
fi

# ============Cleanup============
cprintstart
cprint "Cleaning up ..."

if [[ $(dryrun) == false ]]; then
    if [[ $(system_privileges) == true ]]; then
        mapt-get autoremove
    else
        cprint "Variables have been added to $(project_resource)."
    fi

    cprint "All done! You should now be able to install $(project_name)."
else
    cprint "This was a dryrun, so it didn't install anything."
fi

timer
cd ${GLOBAL_WD}

