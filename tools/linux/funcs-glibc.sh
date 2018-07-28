#!/bin/bash
# 
# Install glibc.
# 

SCRIPT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${SCRIPT_ROOT}/funcs.sh
source ${SCRIPT_ROOT}/funcs-gnu.sh

function glibc_version()
{
    if [[ $(cmdexists ldd) == false ]]; then
        echo 0
    else
        ldd --version | head -1 | awk '{print $5}'
    fi
}

function require_glibc()
{
    cprintstart
    cprint "Verify glibc ${1} ..."

    # Try system installation
    if [[ $(cmdexists ldd) == false ]]; then
        mapt-get install libc6
    fi

    # Check version
    if [[ $(require_new_version $(glibc_version) ${1}) == false ]]; then
        cprint "glibc version $(glibc_version) will be used"
        cprintend
        return
    fi

    # Install from source
    gnu_install_fromsource glibc ${1}

    # Check version
    if [[ $(require_new_version $(glibc_version) ${1}) == false ]]; then
        cprint "glibc version $(glibc_version) will be used"
    else
        if [[ $(cmdexists glibc) == false ]]; then
            cerror "glibc is not installed"
        else
            cerror "glibc version $(glibc_version) will be used but ${1} is required"
        fi
    fi

    cprintend
}


