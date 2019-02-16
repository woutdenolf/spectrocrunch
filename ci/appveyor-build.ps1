# 
# This script will install Appveyor test dependencies.
# 

function appveyor_unittest()
{
    # Build package
    cd $env:APPVEYOR_BUILD_FOLDER
    python setup.py build
    python setup.py sdist
    python setup.py bdist_wheel

    # Install package
    # Pip on Windows: "Requirement already satisfied" when in repo
    $PROJECTNAME = python setup.py name | Select-Object -Last 1
    cd $env:CACHED_FOLDER
    python -m pip install --pre --no-index --find-links="$env:APPVEYOR_BUILD_FOLDER\dist" $PROJECTNAME
}

function appveyor_styletest()
{
    Write-Host "Nothing to build"
}

function main()
{
    if ($env:APPVEYORRUN -eq "unit") {
        appveyor_unittest
    } elseif ($env:APPVEYORRUN -eq "style") {
        appveyor_styletest
    } else {
        Write-Host "Nothing to build"
    }
}

main
