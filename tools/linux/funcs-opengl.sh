#!/bin/bash
# 
# Install pyopengl and drivers.
# 

SCRIPT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${SCRIPT_ROOT}/funcs.sh
source ${SCRIPT_ROOT}/funcs-python.sh


function pyopengl_install()
{
    # Mesa: free OpenGl implementation
    mapt-get install libgl1-mesa-dev
    mapt-get install mesa-common-dev
    # GLU, GLX, GLE: higher level features on top of OpenGl
    mapt-get install libglu1-mesa-dev
    mapt-get install libegl1-mesa
    mapt-get install libgl1-mesa-glx
    pip_install pyopengl
}


function pyopengl_test()
{
    python_get $'try:\n from OpenGL.GL import *\n print("true")\nexcept Exception:\n print("false")'
}


function require_pyopengl()
{
    cprintstart "Require pyopengl"

    # Requirements (for running)
    require_python

    # Check
    if [[ $(pyopengl_test) == true ]]; then
        cprint "Python module \"pyopengl\" is working"
        cprintend "Require pyopengl"
        return
    fi

    # Install
    pyopengl_install

    # Check
    if [[ $(pyopengl_test) == true ]]; then
        cprint "Python module \"pyopengl\" is working"
    else
        cprint "Python module \"pyopengl\" is NOT working"
    fi

    cprintend "Require pyopengl"
}
