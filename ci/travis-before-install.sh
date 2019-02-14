#!/bin/bash
# 
# This script will prepare Travis.
# 

function travis_unittest()
{
    local ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
    mkdir -p $CACHED_FOLDER
    cd $CACHED_FOLDER
    if [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
        source $ROOT/travis-linux-prepare.sh
    fi
}

function main()
{
    if [[ ${TRAVISRUN} == "unit" ]]; then
        travis_unittest
    elif [[ ${TRAVISRUN} == "style" ]]; then
        echo ""
    else
        echo "No tests to be prepared"
    fi
}

main()
