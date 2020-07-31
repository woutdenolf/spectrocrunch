# 
# This script will cleans up Appveyor.
# 

function main()
{
    # Remove virtual environment
    deactivate
    rm $env:VENV_BUILD_DIR -r -fo -ErrorAction Ignore
}

main
