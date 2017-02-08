#!/bin/bash
# 
# Install xraylib on Linux.
# 

# ============Initialize environment============
if [ -z $PYTHONBIN ]; then
    PYTHONBIN=python
fi
if [ -z $PYTHONBINAPT ]; then
    PYTHONBINAPT=python
fi
PYTHON_EXECUTABLE=$(which $PYTHONBIN) # full path
PYTHON_INCLUDE_DIR=`$PYTHONBIN -c "import distutils.sysconfig; print(distutils.sysconfig.get_python_inc());"`
PYTHON_LIBRARY=`$PYTHONBIN -c "import distutils.sysconfig,os; print(os.path.join(distutils.sysconfig.get_config_var('LIBDIR'),distutils.sysconfig.get_config_var('LDLIBRARY')));"`

if [ -z $NOTDRY ]; then
    NOTDRY=true
fi

if [ -z $BUILDSTEP ]; then
    BUILDSTEP=0
fi

if [ -z $TIMELEFT ]; then
    TIMELEFT=true
fi

if [ -z $TIMELIMITED ]; then
    TIMELIMITED=false
fi

if [ -z $SYSTEM_PRIVILIGES ]; then
    if [[ -z "$((sudo -n true) 2>&1)" ]]; then
        export SYSTEM_PRIVILIGES=true 
    else
        export SYSTEM_PRIVILIGES=false
    fi
fi

# ============Install dependencies============
echo -e "${hcol}Install SimpleElastix dependencies ...${ncol}"
if [[ $NOTDRY == true && $SYSTEM_PRIVILIGES == true ]]; then
    sudo -E apt-get -y install swig
    sudo -E apt-get install $PYTHONBINAPT-dev
fi

dpkg --compare-versions "$(cmake --version | head -1 | awk '{print $3}')" "lt" "3.0.0"
if [ $? = 0 ]; then
    RESTORE_WD=$(pwd)
    SCRIPT_WD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
    source $SCRIPT_WD/install-cmake.sh
    cd $RESTORE_WD
else
    BUILDSTEP=$(( $BUILDSTEP+1 ))
fi

# ============Install simpleelastix============
if [ ! -f simpleelastix/build/SimpleITK-build/Wrapping/Python/Packaging/setup.py ]; then
    echo -e "${hcol}Download SimpleElastix ...${ncol}"
    mkdir -p simpleelastix
    cd simpleelastix
    if [[ $NOTDRY == true && ! -d SimpleElastix ]]; then
        git clone https://github.com/kaspermarstal/SimpleElastix SimpleElastix
    fi
    mkdir -p build
    cd build

    if [[ $TIMELEFT == true && ! -f Makefile ]]; then
        echo -e "${hcol}Configure SimpleElastix ...${ncol}"
        if [[ $NOTDRY == true ]]; then
            # http://simpleelastix.readthedocs.io/GettingStarted.html#manually-building-on-linux
            #
            # TODO: USE_SYSTEM_ITK, USE_SYSTEM_ELASTIX
            if [[ $SYSTEM_PRIVILIGES == true ]]; then
                cmake -DBUILD_EXAMPLES:BOOL=OFF \
                            -DBUILD_SHARED_LIBS:BOOL=OFF \
                            -DBUILD_TESTING:BOOL=OFF \
                            -DUSE_SYSTEM_SWIG:BOOL=ON \
                            -DPYTHON_EXECUTABLE:FILEPATH=$PYTHON_EXECUTABLE \
                            -DPYTHON_INCLUDE_DIR:PATH=$PYTHON_INCLUDE_DIR \
                            -DPYTHON_LIBRARY:FILEPATH=$PYTHON_LIBRARY \
                            -DWRAP_CSHARP:BOOL=OFF \
                            -DWRAP_JAVA:BOOL=OFF \
                            -DWRAP_LUA:BOOL=OFF \
                            -DWRAP_PYTHON:BOOL=ON \
                            -DWRAP_R:BOOL=OFF \
                            -DWRAP_RUBY:BOOL=OFF \
                            -DWRAP_TCL:BOOL=OFF \
                            ../SimpleElastix/SuperBuild
            else
                mkdir -p $HOME/.local
                cmake -DCMAKE_INSTALL_PREFIX:PATH="$HOME/.local" \
                            -DBUILD_EXAMPLES:BOOL=OFF \
                            -DBUILD_SHARED_LIBS:BOOL=OFF \
                            -DBUILD_TESTING:BOOL=OFF \
                            -DUSE_SYSTEM_SWIG:BOOL=ON \
                            -DPYTHON_EXECUTABLE:FILEPATH=$PYTHON_EXECUTABLE \
                            -DPYTHON_INCLUDE_DIR:PATH=$PYTHON_INCLUDE_DIR \
                            -DPYTHON_LIBRARY:FILEPATH=$PYTHON_LIBRARY \
                            -DWRAP_CSHARP:BOOL=OFF \
                            -DWRAP_JAVA:BOOL=OFF \
                            -DWRAP_LUA:BOOL=OFF \
                            -DWRAP_PYTHON:BOOL=ON \
                            -DWRAP_R:BOOL=OFF \
                            -DWRAP_RUBY:BOOL=OFF \
                            -DWRAP_TCL:BOOL=OFF \
                            ../SimpleElastix/SuperBuild
            fi

            if [[ $TIMELIMITED == true ]]; then
                    TIMELEFT=false
            fi
        fi
        
        BUILDSTEP=$(( $BUILDSTEP+1 ))
    fi

    if [[ $TIMELEFT == true ]]; then
        echo -e "${hcol}Build SimpleElastix ...${ncol}"
        OMP_NUM_THREADS=2
        if [[ $NOTDRY == true ]]; then
            make -s -j2
            if [[ $TIMELIMITED == true ]]; then
                    TIMELEFT=false
            fi
        fi
        BUILDSTEP=$(( $BUILDSTEP+1 ))
    fi
else
    cd simpleelastix/build
    BUILDSTEP=$(( $BUILDSTEP+2 ))
fi

if [[ $TIMELEFT == true ]]; then
    echo -e "${hcol}Install SimpleElastix ...${ncol}"
    if [[ $NOTDRY == true ]]; then
        cd ./SimpleITK-build/Wrapping/Python/Packaging
        $PYTHONBIN setup.py install
    fi

    BUILDSTEP=$(( $BUILDSTEP+1 ))
fi
