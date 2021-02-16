#!/bin/bash
# 
# Create Travis pre-builds.
# 


# ============Usage============
show_help()
{
  echo "
        Usage: travis-linux-prebuild  -v version [-d] [-p]

        -v version      Python version to be used (2, 3, 2.7, 3.5, ...).
        -d              Dry run.
        -p              Try downloading a pre-build.

        For Example: ./travis-linux-prebuild -v 3 -d -p

        -h              Help
       "
}

# ============Parse script arguments============
OPTIND=0
ARG_DRY=false
ARG_PYTHONV=""
ARG_RET=-1
ARG_PREBUILD=false
while getopts "v:hdp" opt; do
  case $opt in
    h)
      show_help
      ARG_RET=0
      ;;
    d)
      ARG_DRY=true
      ;;
    p)
      ARG_PREBUILD=true
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
    local GLOBAL_SCRIPT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
    source ${GLOBAL_SCRIPT_ROOT}/funcs-travis-linux-prebuild.sh

    dryrun reset ${ARG_DRY}

    if [[ $(travis_check_platform "xenial") == false ]]; then
        echo "Run on the same platfrom as Travis"
        travis_cleanup_python
        return 1
    fi

    if [[ ${ARG_PYTHONV::1} == 3 ]];then
        ARG_PYTHONV=3.5.6
    else
        ARG_PYTHONV=2.7.15
    fi
    travis_init_python ${ARG_PYTHONV}
    if [[ $? != 0 ]]; then
        travis_cleanup_python
        return 1
    fi

    travis_init_cmake 3.12.4
    if [[ $? != 0 ]]; then
        travis_cleanup_python
        return 1
    fi

    if [[ ${ARG_PREBUILD} == false ]]; then
        travis_prepare
        if [[ $? != 0 ]]; then
            travis_cleanup_python
            return 1
        fi
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

if [[ ${ARG_RET} == -1 ]]; then
    if [[ -z ${ARG_PYTHONV} ]]; then
        show_help
        ARG_RET=1
    else
        main
    fi
fi

