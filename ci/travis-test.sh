#!/bin/bash
# 
# This script will prepare Travis.
# 

function travis_unittest()
{
    cd ${HOME}
    #python -m ${PROJECTNAME}.tests.test_all --log=error
    PYOPENCL_COMPILER_OUTPUT=1
    python -m unittest -v ${PROJECTNAME}.tests.test_all.test_suite
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
        echo "No tests to be run"
    fi
}

main
