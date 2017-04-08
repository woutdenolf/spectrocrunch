#!/bin/bash
# 
# Build ITK.
# 

# ============Initialize environment============
SCRIPT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $SCRIPT_ROOT/funcs.sh
initEnv

# ============Build ITK============
if [[ $TIMELEFT == true && ! -f ITK/build/Makefile ]]; then
    cprint "Download ITK ..."
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

    cprint "Configure ITK ..."
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
            mkdir -p $SPECTROCRUNCHLOCAL
            cmake -DCMAKE_INSTALL_PREFIX:PATH="$SPECTROCRUNCHLOCAL" $CMAKE_PARAMS ../ITK
        fi
    fi

    BUILDSTEP=$(( $BUILDSTEP+1 ))

    if [[ $TIMELEFT == true ]]; then
        cprint "Build ITK ..."
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
    cprint "Install ITK ..."
    if [[ $NOTDRY == true ]]; then
        mexec "make install -s"

        if [[ $SYSTEM_PRIVILIGES == false ]]; then
            addProfile $SPECTROCRUNCHRC "# Installed ITK: $SPECTROCRUNCHLOCALSTR"
            addLibPath $SPECTROCRUNCHLOCAL/lib
            addLibPathProfile $SPECTROCRUNCHRC "$SPECTROCRUNCHLOCALSTR/lib"
        fi

        # Variable needed by simpleelastix:
        ITK_DIR=$SPECTROCRUNCHLOCAL/lib/cmake
        ITK_DIR=$ITK_DIR/$(find $ITK_DIR -maxdepth 1 -type d -name 'ITK-*' -printf %f -quit)
        cprint "ITK directory: $ITK_DIR"
    fi

    BUILDSTEP=$(( $BUILDSTEP+1 ))
fi

BUILDSTEPS=$(( $BUILDSTEPS+1 ))
cd $INSTALL_WD


