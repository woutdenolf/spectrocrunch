# 
# This script will install Appveyor test dependencies.
# 

function appveyor_unittest()
{
    # Install build dependencies
    cd $env:CACHED_FOLDER
    if ($env:OS -eq "Windows_NT") {
        & "$env:APPVEYOR_BUILD_FOLDER\tools\windows-install-deps.ps1" -y -u -s
    }
}

function appveyor_styletest()
{
    # Install build dependencies
    cd $env:CACHED_FOLDER
    python -m pip install flake8
}

function main()
{
    if ($env:APPVEYORRUN -eq "unit") {
        appveyor_unittest
    } elseif ($env:APPVEYORRUN -eq "style") {
        appveyor_styletest
    } else {
        Write-Host "No dependencies to be installed"
    }

    # Print Python info
    Write-Host "Python information:"
    python $env:APPVEYOR_BUILD_FOLDER\ci\info_platform.py
    python -m pip list
}

main
