# 
# This script will prepare Appveyor.
# 


function main()
{
    # Environment
    [Environment]::SetEnvironmentVariable("Path", "$env:PYTHON_DIR;$env:PYTHON_DIR\Scripts;$env:Path", "Process")
    Write-Host "Build worker environment variables:" -ForegroundColor Magenta
    Get-ChildItem env:

    # Python virtual environment
    $env:PYTHONWARNINGS="ignore"
    python -m ensurepip
    python -m pip install --upgrade pip
    python -m pip install --upgrade setuptools
    python -m pip install --upgrade virtualenv
    virtualenv --clear $env:VENV_BUILD_DIR
    & "$env:VENV_BUILD_DIR\Scripts\activate.ps1"

    # Create folder for custom cache
    mkdir $env:CACHED_FOLDER -Force
    cd $env:CACHED_FOLDER
}

main
