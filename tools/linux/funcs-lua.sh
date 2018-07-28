#!/bin/bash
# 
# Install lua.
# 

SCRIPT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${SCRIPT_ROOT}/funcs.sh

function lua_version()
{
    if [[ $(cmdexists lua) == false ]]; then
        echo 0
    else
        lua -v 2>&1 | awk '{print $2}'
    fi
}

