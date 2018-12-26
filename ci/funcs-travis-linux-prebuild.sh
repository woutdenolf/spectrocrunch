#!/bin/bash
# 
# Functions to create Travis pre-builds.
# 

GLOBAL_SCRIPT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${GLOBAL_SCRIPT_ROOT}/../tools/linux/funcs.sh
source ${GLOBAL_SCRIPT_ROOT}/../tools/linux/funcs-python.sh

function travis_check_platform()
{
    local platform=$(lsb_release -c | awk '{print $2}')
    if [[ ${platform} == ${1} ]]; then
        echo true
    else
        echo false
    fi
}

function travis_build_folder()
{
    echo "/home/travis"
}


function travis_cached_folder()
{
    local cached_folder=$(travis_build_folder)/cached
    mkdir -p ${cached_folder}
    echo "${cached_folder}"
}


function travis_pybuild_folder()
{
    local pybuild_folder=$(travis_cached_folder)/$(project_name)_$(python_full_version)
    mkdir -p ${pybuild_folder}
    echo "${pybuild_folder}"
}


function travis_depbuild_folder()
{
    local pybuild_folder=$(travis_cached_folder)/$(python_depdir)
    mkdir -p ${pybuild_folder}
    echo "${pybuild_folder}"
}


function travis_venv()
{
    echo $(travis_build_folder)/virtualenv/python$(python_full_version)
}


function travis_init_python()
{
    local restorewd=$(pwd)
    cd $(travis_cached_folder)

    local pythonv=${1}

    if [[ $(system_privileges) == false ]]; then
        sudo -s "exit"
    fi

    if [[ $(dryrun) == false ]]; then
        # Require python version
        install_systemwide reset true
        python_virtualenv_deactivate
        require_python ${pythonv}
        if [[ $? != 0 ]]; then
            cd ${restorewd}
            return 1
        fi

        # Require pip
        install_systemwide reset false
        require_pip
        pip_install virtualenv

        # Activate virtual environment
        virtualenv $(travis_venv)
        source $(travis_venv)/bin/activate
    fi

    cd ${restorewd}
}


function travis_install_dependencies()
{
    local restorewd=$(pwd)
    cd $(travis_cached_folder)

    local pythonv=${1}

    if [[ $(dryrun) == true ]]; then
        source $(project_folder)/tools/linux-install-deps.sh -v ${pythonv} -u -d -x
    else
        source $(project_folder)/tools/linux-install-deps.sh -v ${pythonv} -u -y -x
    fi

    local ret=$?
    cd ${restorewd}
    return ${ret}
}


function travis_prepare()
{
    local restorewd=$(pwd)
    cd $(travis_cached_folder)

    if [[ $(dryrun) == false ]]; then
        source $(project_folder)/ci/travis-linux-prepare.sh
    fi

    cd ${restorewd}
}


function travis_build_project()
{
    local restorewd=$(pwd)
    cd $(project_folder)

    cprintstart
    cprint "Build project ..."

    if [[ $(dryrun) == false ]]; then
        $(python_bin) setup.py build_py --build-lib=$(travis_pybuild_folder)/build/lib -f
        $(python_bin) setup.py build_doc --build-dir=$(travis_pybuild_folder)/build/doc
        $(python_bin) setup.py sdist --dist-dir=$(travis_pybuild_folder)/dist
        $(python_bin) setup.py bdist_wheel --dist-dir=$(travis_pybuild_folder)/dist
    fi

    cd ${restorewd}
}


function travis_test_project()
{
    local restorewd=$(pwd)
    cd $(travis_build_folder)

    cprintstart
    cprint "Install/test $(project_name) ..."

    if [[ $(dryrun) == false ]]; then
        cprint "Uninstall $(project_name) ..."
        pip_uninstall $(project_name)

        cprint "Install $(project_name) ..."
        pip_install --pre --no-index --no-cache --find-links=$(travis_pybuild_folder)/dist/ $(project_name)

        cprint "Test $(project_name) ..."
        $(python_bin) -m $(project_name).tests.test_all
    fi

    cd ${restorewd}
}


function travis_cleanup_python()
{
    if [[ $(dryrun) == false ]]; then
        python_virtualenv_deactivate
    fi
}


function travis_pack_prebuild()
{
    local restorewd=$(pwd)
    cd $(travis_cached_folder)

    local filename=$(project_name).travis.python$(python_full_version).tgz
    cprintstart
    cprint "Pack ${filename} ..."

    if [[ $(dryrun) == false ]]; then
        rm -f ${filename}
        tar -czf ${filename} $(basename $(travis_depbuild_folder))

        #require_web_essentials
        #local link=$(curl --upload-file ./${filename} https://transfer.sh/${filename})
        #echo "Link for download: ${link}"
    fi

    cd ${restorewd}
}
