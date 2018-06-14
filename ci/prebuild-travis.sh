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
while getopts "v:uyhd" opt; do
  case $opt in
    h)
      show_help
      return 0
      ;;
    d)
      ARG_DRY=true
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

PROJECT="$(basename  $( cd "$( dirname "${BASH_SOURCE[0]}" )"/.. && pwd ))"
if [[ -z ${ARG_PYTHONV} ]]; then
    show_help
    return 1
fi

# ============Initialize environment============
BUILD_FOLDER=/home/travis
CACHED_FOLDER=${BUILD_FOLDER}/cached
PROJECT_FOLDER=${BUILD_FOLDER}/${PROJECT}

mkdir -p ${CACHED_FOLDER}

cd ${BUILD_FOLDER}
source ${PROJECT_FOLDER}/tools/linux/funcs-python.sh
$(dryrun reset ${ARG_DRY})
$(install_systemwide reset false)

if [[ $(system_privileges) == false ]]; then
    sudo -s
fi

if [[ $(dryrun) == false ]]; then
    if [[ $(cmdexists deactivate) == true ]]; then
        deactivate
    fi
    require_python ${PYTHONV}
    require_pip
    return 0

    pip_install virtualenv
    virtualenv python${PYTHONV}
    source python${PYTHONV}/bin/activate
fi

dryrun
return 0

# ============Install dependencies============
PREBUILD_FOLDER=${CACHED_FOLDER}/$(python_full_version)
mkdir -p ${PREBUILD_FOLDER}
cd ${PREBUILD_FOLDER}
if [[ $(dryrun) == true ]]; then
    . ${PROJECT_FOLDER}/tools/prepare_install-linux.sh -v ${ARG_PYTHONV} -u -d
else
    . ${PROJECT_FOLDER}/tools/prepare_install-linux.sh -v ${ARG_PYTHONV} -u
fi
if [[ $? != 0 ]]; then
    return 1
fi

# ============Build project============
PYBUILD_FOLDER=${CACHED_FOLDER}/$(python_full_version)/${PROJECT}
mkdir -p ${PYBUILD_FOLDER}
cd ${PYBUILD_FOLDER}
if [[ $(dryrun) == false ]]; then
    python_bin ${PROJECT_FOLDER}/setup.py build

    mapt-get pandoc # nbsphinx
    pip_install -r ${PROJECT_FOLDER}/requirements-dev.txt
    python_bin ${PROJECT_FOLDER}/setup.py build_doc

    python_bin ${PROJECT_FOLDER}/setup.py sdist bdist_wheel
fi

# ============Test build============
cd ${BUILD_FOLDER}
if [[ $(dryrun) == false ]]; then
    pip_install --pre --no-index --find-links=${PROJECT_FOLDER}/dist/ ${PROJECT}
    python_bin -m ${PROJECT}.tests.test_all
    deactivate
fi

# ============Pack result============
#tar -czvf spectrocrunch.travis.python${PYTHONV}.tgz simpleelastix

