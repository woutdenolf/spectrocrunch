function appveyor_download_prebuild()
{
    # Download pre-build libraries
    $local:pythonv = python -c "import sys;t='{v[0]}.{v[1]}.{v[2]}'.format(v=list(sys.version_info[:3]));print(t)"

    $local:dep_folder = "dep_$local:pythonv"
}

function main()
{
    appveyor_download_prebuild
}

main