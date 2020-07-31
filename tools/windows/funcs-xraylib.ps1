# 
# Install xraylib on Windows.
# 

. $PSScriptRoot\funcs-python.ps1


function xraylib_url()
{
    return "http://lvserver.ugent.be/xraylib/"
}


function xraylib_latest($extension) {
    $local:pattern = version_pattern $null 3
    $local:pattern += $extension
    $local:mversion = (0,0,0)
    $local:dllink = $null
    $local:url = $(xraylib_url)
    $local:content = Invoke-WebRequest $(xraylib_url)
    foreach ($link in $local:content.Links) {
        $local:m = [regex]::match($link.href,$local:pattern)
        if ($local:m.Success) {
            $local:tmp = version_intarray $m.Groups[1].Value 3
            if ((version_intarray_higher $local:mversion $local:tmp)) {
                $local:mversion = $local:tmp
                $local:dllink = "$local:url/"+$link.href
            }
        }
    }
    return $local:dllink
}


function xraylib_install()
{
    $local:restorewd = (Get-Location).Path
    mkdir xraylib -Force
    cd xraylib

    if ((python_arch) -eq 64) {
        $local:extension = "-win64.exe"
        $local:affix = ""
        $local:name = "xraylib (64-bit)"
    } else {
        $local:extension = "-win32.exe"
        $local:affix = "-32"
        $local:name = "xraylib (32-bit)"
    }

    cprint "Download xraylib ..."
    if (!(dryrun)) {
        require_web_essentials
        $local:filename = @(gci "*$local:extension")[0]
        if ($local:filename -eq $null) {
            $local:url = xraylib_latest $local:extension
            cprint $local:url
            $local:filename = download_file $local:url ""
        }
    }

    $local:path = joinPath (project_prefix) "xraylib$affix"
    cprint "Installing $name in $path ..."
    if (!(dryrun) -and $local:filename -ne $null) {
        $arguments = @()
        $arguments += "/silent"
        $arguments += "/dir=""$path"""
        install_exe $filename $arguments
        prependBinPath (joinPath $path "bin")
        if (cmdexists xraylib) {
            cprint "xraylib is installed."
        } else {
            cerror "xraylib is not installed."
        }
    }

    cd $local:restorewd
}


function require_xraylib()
{
    cprintstart
    cprint "Verify xraylib ..."

    # Requirements (for running)
    require_python

    # Check
    if (python_hasmodule "xraylib") {
        cprint "Python module ""xraylib"" is installed"
        cprintend
        return
    }

    # Install xraylib
    xraylib_install

    # Check
    if (python_hasmodule "xraylib") {
        cprint "Python module ""xraylib"" is installed"
    } else {
        cprint "Python module ""xraylib"" is NOT installed"
    }

    cprintend
}