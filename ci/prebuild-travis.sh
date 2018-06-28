#!/bin/bash
# 
# Create Travis pre-builds.
# 


# ============Usage============
show_help()
{
  echo "
        Usage: prebuild_travis  -v version [-d]

        -v version      Python version to be used (2, 3, 2.7, 3.5, ...).
        -d              Dry run.

        For Example: ./prebuild_travis -v 3 -d

        -h              Help
       "
}

# ============Parse script arguments============
OPTIND=0
ARG_DRY=false
ARG_PYTHONV=""
ARG_RET=-1
while getopts "v:hd" opt; do
  case $opt in
    h)
      show_help
      ARG_RET=0
      ;;
    d)
      ARG_DRY=true
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

if [[ -z ${ARG_PYTHONV} ]]; then
    show_help
    ARG_RET=1
fi


function main()
{
    if [[ ${ARG_RET} != -1 ]]; then
        return ${ARG_RET}
    fi

    # ============Initialize environment============
    local GLOBAL_SCRIPT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
    source ${GLOBAL_SCRIPT_ROOT}/funcs-prebuild-travis.sh
    dryrun reset ${ARG_DRY}

    travis_init_python ${ARG_PYTHONV}
    if [[ $? != 0 ]]; then
        return 1
    fi

    # ============Dependencies============
    travis_install_dependencies ${ARG_PYTHONV}
    if [[ $? != 0 ]]; then
        travis_cleanup_python
        return 1
    fi

    # ============Build============
    travis_build_project

    # ============Test============
    travis_test_project

    # ============Pack============
    travis_pack_prebuild

    # ============Cleanup============
    travis_cleanup_python
}

main

