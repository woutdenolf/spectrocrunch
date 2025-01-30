#!/bin/bash
# 
# This script will run Travis tests.
# 

function travis_unittest()
{
    cd ${TRAVIS_BUILD_DIR}
    PROJECTNAME=`python setup.py name|tail -1`

    cd ${HOME}
    PYOPENCL_COMPILER_OUTPUT=1
    #python -m ${PROJECTNAME}.tests.test_all --log=error
    python -m unittest -v ${PROJECTNAME}.tests.test_all.main_test_suite
}

function travis_styletest()
{
    cd ${TRAVIS_BUILD_DIR}
    python -m flake8 spectrocrunch
}

function main()
{
    if [[ ${TRAVISRUN} == "unit" ]]; then
        travis_unittest
    elif [[ ${TRAVISRUN} == "style" ]]; then
        travis_styletest
    else
        echo "No tests to run"
    fi
}

main
