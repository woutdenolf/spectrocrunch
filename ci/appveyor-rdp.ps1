# 
# This script will block Appveyor and allow RDP.
# 

function main()
{
    # Keep environment for RDP login
    . "$env:APPVEYOR_BUILD_FOLDER\tools\windows\funcs-install.ps1"
    copyEnv "Process" "User"

    # Blocking RDP start
    $blockRdp = $true
    iex ((new-object net.webclient).DownloadString('https://raw.githubusercontent.com/appveyor/ci/master/scripts/enable-rdp.ps1'))
}

main
