#
# Basic powershell requirements
#

# ============cprint============
# Description: output to stdout in color
# Usage: cprint "..."
function cprint() 
{
    Write-Host -ForegroundColor Magenta $args
}


# ============cerror============
# Description: output to stdout in color
# Usage: cerror "..."
function cerror()
{
    Write-Host -ForegroundColor Red $args
}


# Check version
if ($PSVersionTable.PSVersion.Major -lt 3) {
    cerror "To run this script you need to download WMF 3.0 or higher:"
    cerror " https://docs.microsoft.com/en-us/powershell/wmf"
    $global:ErrorActionPreference = "Stop"
    throw "Download WMF 3.0 or higher"
}
