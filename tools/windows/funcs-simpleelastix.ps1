# 
# Install simpleelastix on Windows.
# 

. $PSScriptRoot\funcs-python.ps1
. $PSScriptRoot\funcs-git.ps1
. $PSScriptRoot\funcs-cmake.ps1

function simpleelastix_build_dependencies()
{
    require_cmake 3
    python_init_compiler
}


function simpleelastix_install_fromsource()
{
    if (!(path_exists simpleelastix) -and $global:ARG_SKIPLONG) {
        cprint "Skipping simpleelastix installation"
        return
    }

    $local:restorewd = (Get-Location).Path

    cprint "Download SimpleElastix ..."
    mkdir simpleelastix -Force
    cd simpleelastix
    if (!(path_exists simpleelastix) -and !(dryrun)) {
        require_web_essentials
        require_git
        git clone https://github.com/kaspermarstal/SimpleElastix SimpleElastix
    }
    mkdir build -Force
    cd build

    $local:prefix = joinPath $(project_prefix) "simpleelastics"
    $local:outdir = joinPath (Get-Location).path "SimpleITK-build\Wrapping\Python\Packaging"
    if (!(dryrun) -and !(file_exists $(joinPath $local:outdir "setup.py")))
    {
        cprint "Configure SimpleElastix ..."
        if (!(dryrun)) {
            simpleelastix_build_dependencies

            # https://simpleelastix.readthedocs.io/GettingStarted.html#compiling-on-windows
            $local:CMAKE_PARAMS=@()
            if ((cmdexists cl) -and (current_msc_ver) -eq 1500) {
                $local:_COMPILER = (Get-Command cl).definition
                $local:CMAKE_PARAMS += "-DCMAKE_C_COMPILER:FILEPATH=""$local:_COMPILER"""
                $local:CMAKE_PARAMS += "-DCMAKE_CXX_COMPILER:FILEPATH=""$local:_COMPILER"""
            } else {
                $local:_COMPILER = msc_cmake_info
                $local:_G = $local:_COMPILER["generator"]
                $local:_A = $local:_COMPILER["arch"]
                $local:_T = $local:_COMPILER["toolset"]
                if ($local:_G -ne $null) {
                    $local:CMAKE_PARAMS += "-G ""$local:_G"""
                }
                if ($local:_A -ne $null) {
                    $local:CMAKE_PARAMS += "-A $local:_A"
                }
                if ($local:_T -ne $null) {
                    $local:CMAKE_PARAMS += "-T $local:_T"
                }
            }
            $local:CMAKE_PARAMS += "-DCMAKE_INSTALL_PREFIX:PATH=""$local:prefix"""
            $local:CMAKE_PARAMS += "-DBUILD_EXAMPLES:BOOL=OFF"
            $local:CMAKE_PARAMS += "-DBUILD_SHARED_LIBS:BOOL=OFF"
            $local:CMAKE_PARAMS += "-DBUILD_TESTING:BOOL=OFF"
            $local:CMAKE_PARAMS += "-DSimpleITK_USE_SYSTEM_SWIG:BOOL=OFF"
            $local:CMAKE_PARAMS += "-DSimpleITK_USE_SYSTEM_LUA:BOOL=OFF"
            $local:CMAKE_PARAMS += "-DSimpleITK_USE_SYSTEM_VIRTUALENV:BOOL=OFF"
            $local:CMAKE_PARAMS += "-DSimpleITK_USE_SYSTEM_ELASTIX:BOOL=OFF"
            $local:CMAKE_PARAMS += "-DSimpleITK_USE_SYSTEM_ITK:BOOL=OFF"
            $local:CMAKE_PARAMS += "-DPYTHON_EXECUTABLE:FILEPATH=$(python_bin)"
            $local:CMAKE_PARAMS += "-DPYTHON_INCLUDE_DIR:PATH=$(python_include)"
            $local:CMAKE_PARAMS += "-DPYTHON_LIBRARY:FILEPATH=$(python_lib)"
            $local:CMAKE_PARAMS += "-DWRAP_DEFAULT:BOOL=OFF"
            $local:CMAKE_PARAMS += "-DWRAP_PYTHON:BOOL=ON"
            $local:CMAKE_PARAMS += "..\SimpleElastix\SuperBuild"

            mkdir $local:prefix -Force
            Start-Process "cmake" -ArgumentList $(@("-LAH")+$local:CMAKE_PARAMS) -PassThru -NoNewWindow -Wait
            #cmake $local:CMAKE_PARAMS
        }

        cprint "Build SimpleElastix ..."
        if (!(dryrun)) {
            echo msbuild /p:configuration=release ALL_BUILD.vcxproj
        }
    }

    cprint "Install SimpleElastix ..."
    if (!(dryrun)) {
        cd $local:outdir
        #python_exec setup.py install

        #addProfile $(project_resource) "# Installed simpleelastix: ${prefixstr}"
        #addLibPath ${prefix}/lib
        #addLibPathProfile $(project_resource) "${prefixstr}/lib"
    }

    cd $local:restorewd
}


function require_simpleelastix()
{
    cprintstart
    cprint "Verify simpleelastix ..."

    # Requirements (for running)
    require_python

    # Check
    if (python_hasmodule "SimpleITK") {
        cprint "Python module ""simpleelastix"" is installed"
        cprintend
        return
    }

    # Install from source
    simpleelastix_install_fromsource

    # Check
    if (python_hasmodule "SimpleITK") {
        cprint "Python module ""simpleelastix"" is installed"
    } else {
        cprint "Python module ""simpleelastix"" is NOT installed"
    }

    cprintend
}