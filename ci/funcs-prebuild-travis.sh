#!/bin/bash
# 
# Functions to create Travis pre-builds.
# 

GLOBAL_SCRIPT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${GLOBAL_SCRIPT_ROOT}/../tools/linux/funcs.sh
source ${GLOBAL_SCRIPT_ROOT}/../tools/linux/funcs-python.sh

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


function travis_prebuild_folder()
{
    local prebuild_folder=$(travis_cached_folder)/build_$(python_full_version)
    mkdir -p ${prebuild_folder}
    echo "${prebuild_folder}"
}


function travis_pybuild_folder()
{
    local pybuild_folder=$(travis_prebuild_folder)/$(project_name)
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
        # Install python version
        install_systemwide reset true
        python_virtualenv_deactivate
        require_python ${pythonv}
        if [[ $? != 0 ]]; then
            cd ${restorewd}
            return 1
        fi

        # Activate virtual environment
        install_systemwide reset false
        require_pip
        pip_install virtualenv
        virtualenv $(travis_venv)
        source $(travis_venv)/bin/activate
    fi

    cd ${restorewd}
}


function travis_install_dependencies()
{
    local restorewd=$(pwd)
    cd $(travis_prebuild_folder)

    local pythonv=${1}

    if [[ $(dryrun) == true ]]; then
        . $(project_folder)/tools/prepare_install-linux.sh -v ${pythonv} -u -d -x
    else
        . $(project_folder)/tools/prepare_install-linux.sh -v ${pythonv} -u -x
    fi

    local ret=$?
    cd ${restorewd}
    return ${ret}
}


function travis_build_project()
{
    local restorewd=$(pwd)
    cd $(travis_pybuild_folder)

    if [[ $(dryrun) == false ]]; then
        $(python_bin) $(project_folder)/setup.py build

        $(python_bin) $(project_folder)/setup.py build_doc

        $(python_bin) $(project_folder)/setup.py sdist bdist_wheel
    fi

    cd ${restorewd}
}


function travis_test_project()
{
    local restorewd=$(pwd)
    cd $(travis_build_folder)

    if [[ $(dryrun) == false ]]; then
        pip_install --pre --no-index --find-links=$(project_folder)/dist/ $(project_name)
        $(python_bin) -m $(project_name).tests.test_all
        deactivate
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
    cd $(travis_build_folder)

    if [[ $(dryrun) == false ]]; then
        local filename=$(project_name).travis.python$(python_version).tgz
        #tar -czvf ${file} simpleelastix

        #require_web_essentials
        #local link=$(curl --upload-file ./${file} https://transfer.sh/${file})
        #echo "Link for download: ${link}"
    fi

    cd ${restorewd}
}




