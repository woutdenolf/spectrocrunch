#!/bin/bash
# 
# Build ITK.
# 

# ============Initialize environment============
if [ -z $NOTDRY ]; then
    NOTDRY=true
fi

if [ -z $BUILDSTEP ]; then
    BUILDSTEP=0
    BUILDSTEPS=0
fi

if [ -z $SYSTEM_PRIVILIGES ]; then
    if [[ -z "$((sudo -n true) 2>&1)" ]]; then
        export SYSTEM_PRIVILIGES=true 
    else
        export SYSTEM_PRIVILIGES=false
    fi
fi

# ============Build ITK============
if [[ $TIMELEFT == true && ! -f ITK/build/Makefile ]]; then
    echo -e "${hcol}Download ITK ...${ncol}"
    mkdir -p ITK
    cd ITK
    if [[ $NOTDRY == true && ! -d ITK ]]; then
        git clone https://itk.org/ITK.git ITK
        cd ITK
        # TODO: fixed version used by the simpleelastix superbuild
        git reset --hard 80178ae516db87aa50d2883c4e090da92c4a502d
        cd ..
    fi
    mkdir -p build
    cd build

    echo -e "${hcol}Configure ITK ...${ncol}"
    if [[ $NOTDRY == true ]]; then
        CMAKE_PARAMS="-DBUILD_SHARED_LIBS=OFF \
                      -DModule_ITKReview=ON \
                      -DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=ON
                      -DCMAKE_VISIBILITY_INLINES_HIDDEN:BOOL=ON
                      -DBUILD_EXAMPLES:BOOL=OFF
                      -DBUILD_TESTING:BOOL=OFF
                      -DBUILD_SHARED_LIBS:BOOL=OFF
                      -DCMAKE_SKIP_RPATH:BOOL=ON
                      -DITK_LEGACY_REMOVE:BOOL=ON
                      -DITK_BUILD_DEFAULT_MODULES:BOOL=ON
                      -DModule_ITKReview:BOOL=ON
                      -DITK_USE_GIT_PROTOCOL:BOOL=OFF
                      -DITK_WRAP_float:BOOL=ON
                      -DITK_WRAP_unsigned_char:BOOL=ON
                      -DITK_WRAP_signed_short:BOOL=ON
                      -DITK_WRAP_unsigned_short:BOOL=ON
                      -DITK_WRAP_complex_float:BOOL=ON
                      -DITK_WRAP_vector_float:BOOL=ON
                      -DITK_WRAP_covariant_vector_float:BOOL=ON
                      -DITK_WRAP_rgb_signed_short:BOOL=ON
                      -DITK_WRAP_rgb_unsigned_char:BOOL=ON
                      -DITK_WRAP_rgb_unsigned_short:BOOL=ON
                      -DITK_WRAP_PYTHON:BOOL=OFF
                      -DINSTALL_WRAP_ITK_COMPATIBILITY:BOOL=OFF"
        if [[ $SYSTEM_PRIVILIGES == true ]]; then
            cmake $CMAKE_PARAMS ../ITK
        else
            mkdir -p $HOME/.local
            cmake -DCMAKE_INSTALL_PREFIX:PATH="$HOME/.local" $CMAKE_PARAMS ../ITK
        fi
    fi

    BUILDSTEP=$(( $BUILDSTEP+1 ))

    if [[ $TIMELEFT == true ]]; then
        echo -e "${hcol}Build ITK ...${ncol}"
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
    cd ITK/build
    BUILDSTEP=$(( $BUILDSTEP+2 ))
fi

BUILDSTEPS=$(( $BUILDSTEPS+2 ))

if [[ $TIMELEFT == true ]]; then
    echo -e "${hcol}Install ITK ...${ncol}"
    if [[ $NOTDRY == true ]]; then
        if [[ $SYSTEM_PRIVILIGES == true ]]; then
            sudo -E make install -s
            export ITK_DIR=/usr/local/lib/cmake
        else
            make install -s
            export ITK_DIR=$HOME/.local/lib/cmake
        fi

        ITK_DIR=$ITK_DIR/$(find $ITK_DIR -maxdepth 1 -type d -name 'ITK-*' -printf %f -quit)

        echo -e "${hcol}ITK directory: $ITK_DIR${ncol}"
    fi

    BUILDSTEP=$(( $BUILDSTEP+1 ))
fi

BUILDSTEPS=$(( $BUILDSTEPS+1 ))



