#!/bin/bash
# 
# This script will install Travis test dependencies.
# 

function travis_unittest()
{
    # Install general build dependencies
    cd $CACHED_FOLDER
    if [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
        source $TRAVIS_BUILD_DIR/tools/linux-install-deps.sh -y -u -x -s
    else
        # Does not install non-pypi libraries: 
        # corresponding tests will be skipped
        type python
        python -m pip install --upgrade pip --user
        python --version
        pip --version
        python -m pip install --upgrade setuptools --user
        python -m pip install --upgrade wheel --user
        python -m pip install pybind11 --user
        python -m pip install -r $TRAVIS_BUILD_DIR/requirements.txt --user
        python -m pip install -r $TRAVIS_BUILD_DIR/requirements-dev.txt --user
    fi
}

function travis_styletest()
{
    cd $CACHED_FOLDER
    python -m pip install flake8
}

function main()
{
    if [[ ${TRAVISRUN} == "unit" ]]; then
        travis_unittest
    elif [[ ${TRAVISRUN} == "style" ]]; then
        travis_styletest
    else
        echo "No dependencies to be installed"
    fi

    # Print Python info
    echo "Python information:"
    python $TRAVIS_BUILD_DIR/ci/info_platform.py
    python -m pip list
}

main
