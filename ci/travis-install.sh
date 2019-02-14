#!/bin/bash
# 
# This script will prepare Travis.
# 

function travis_unittest()
{
    # Install general build dependencies
    cd $CACHED_FOLDER
    if [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
        source $TRAVIS_BUILD_DIR/tools/linux-install-deps.sh -y -u -x -s
    else
        # Does not install non-pypi libraries
        # tests involved will be skipped
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

    # Print Python info
    python $TRAVIS_BUILD_DIR/ci/info_platform.py
    python -m pip list

    # Build package
    cd $TRAVIS_BUILD_DIR
    python setup.py build
    python setup.py sdist bdist_wheel
    if [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
        python setup.py build_doc
    fi

    # Install package
    PROJECTNAME=`python setup.py name|tail -1`
    pip install --pre --no-index --find-links=dist/ ${PROJECTNAME}
}

function travis_styletest()
{
    cd ${TRAVIS_BUILD_DIR}
    flake8 spectrocrunch
}

function main()
{
    if [[ ${TRAVISRUN} == "unit" ]]; then
        travis_unittest
    elif [[ ${TRAVISRUN} == "style" ]]; then
        travis_styletest
    else
        echo "No tests to be prepared"
    fi
}

main()
