#!/bin/bash
# 
# This script will install Travis test dependencies.
# 

function travis_unittest()
{
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
    echo "Nothing to build"
}

function main()
{
    if [[ ${TRAVISRUN} == "unit" ]]; then
        travis_unittest
    elif [[ ${TRAVISRUN} == "style" ]]; then
        travis_styletest
    else
        echo "Nothing to build"
    fi
}

main
