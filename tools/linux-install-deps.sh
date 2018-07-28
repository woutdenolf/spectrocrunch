#!/bin/bash
# 
# This script will install all system and python dependencies.
# 

# ============Usage============
show_help()
{
  echo "
        Usage: linux-install-deps  [-v version] [-y] [-d] [-x]

        -v version      Python version to be used (2, 3, 2.7, 3.5, ...).
        -y              Answer yes to everything.
        -d              Dry run.
        -u              Install for user only.
        -x              Dev install.
        -t              Do not modify bashrc.

        For Example: ./linux-install-deps.sh -v 3 -d

        -h              Help
        
       Environment variables(optional):
       
        PROJECT_PREFIX: directory to install dependencies
        PROJECT_RESOURCE_DIR: directory of the project resource file
       "
}

# ============Parse script arguments============
OPTIND=0
ARG_SYSTEMWIDE=true
ARG_DRY=false
ARG_BASHRC=false
ARG_FORCECHOICE=false
ARG_PYTHONV=""
ARG_DEV=false
ARG_RET=-1
while getopts "v:uyhdx" opt; do
  case $opt in
    h)
      show_help
      ARG_RET=0
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
    y)
      ARG_BASHRC=true
      ;;
    x)
      ARG_DEV=true
      ;;
    v)
      ARG_PYTHONV=${OPTARG}
      ;;
    \?)
      echo "Invalid option: -${OPTARG}. Use -h flag for help." >&2
      ARG_RET=1
      ;;
  esac
done

function main()
{
    if [[ ${ARG_RET} != -1 ]]; then
        return ${ARG_RET}
    fi

    # ============Initialize environment============
    local GLOBAL_WD=$(pwd)
    local GLOBAL_SCRIPT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
    source $GLOBAL_SCRIPT_ROOT/linux/funcs-python.sh
    dryrun reset ${ARG_DRY}
    modify_bashrc reset ${ARG_BASHRC}
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
    install_info

    local ARG_CHOICE
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

    mkdir -p $(python_depdir)
    cd $(python_depdir)

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
}

main

