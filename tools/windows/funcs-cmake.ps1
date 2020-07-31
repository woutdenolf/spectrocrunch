#
# Cmake
#

. $PSScriptRoot\funcs.ps1
. $PSScriptRoot\funcs-install.ps1


function cmake_url()
{
    return "https://cmake.org/files/"
}


function cmake_download_link([string]$version)
{
    $local:pattern1 = version_pattern $version 2
    $local:pattern2 = version_pattern $version 3
    if ((install_arch) -eq 64) {
        $local:lversion,$local:link = get_version_link $(cmake_url) "v$local:pattern1" "$local:pattern2-(?!rc).*x64\.[exmsi]{3}$"
        if ($local:link -ne $null) {
            return $local:lversion,$local:link
        } 
    }
    return get_version_link $(cmake_url) "v$local:pattern1" "$local:pattern2-(?!rc).*\.[exmsi]{3}$" 
}


function cmake_install_fromsource([string]$version)
{
    $local:restorewd = (Get-Location).Path

    cprint "Download cmake ..."
    mkdir cmake -Force
    cd cmake

    $local:lversion,$local:filename = get_local_version $version
    if ($local:filename -eq $null){
        require_web_essentials
        $local:version,$local:link = cmake_download_link $version
        if ($local:link -eq $null) {
            cerror "No cmake $version found for this platform"
        } else {
            if (!(dryrun)) {
                $local:filename = download_file $local:link
            }
        }
    } else {
        $local:version = $local:lversion
    }

    cprint "Installing cmake $local:version ..."
    if (!(dryrun) -and $local:filename -ne $null) {
        $arguments = @{}
        $arguments["msi"] = @()
        $arguments["exe"] = @()

        $arguments["msi"] += "/passive"
        $arguments["exe"] += "/S"

        $local:systemwide = [int]$(install_systemwide)
        $arguments["msi"] += "ALLUSERS=""$local:systemwide"""
        $arguments["exe"] += "InstallAllUsers=$local:systemwide"

        install_any $local:filename $arguments

        $local:installpath = registry-value "HKEY_LOCAL_MACHINE\Software\Kitware\CMake $version" "(Default)"
        if ($local:installpath -ne $null) {
            prependBinPath $(joinPath $local:installpath "bin")
        }
        $local:installpath = registry-value "Registry::HKEY_LOCAL_MACHINE\SOFTWARE\Kitware\CMake" "InstallDir"
        if ($local:installpath -ne $null) {
            prependBinPath $(joinPath $local:installpath "bin")
        }

        if (cmdexists cmake) {
            cprint "cmake is installed."
        } else {
            cerror "cmake is not installed."
        }
    }

    cd $local:restorewd
}

function cmake_version()
{
    if (cmdexists cmake) {
        $local:m = [regex]::match((cmake --version)[0],"([\d\.]+)")
        if ($local:m.Success) {
            return $m.Groups[1].Value
        }
    }
    return $null
}

function require_cmake([AllowNull()][string]$version)
{
    cprintstart
    cprint "Verify cmake $version ..."

    # Check version
    if (!(require_new_version $(cmake_version) $version)) {
        cprint "cmake version $(cmake_version) will be used"
        cprintend
        return
    }

    # Install from source
    cmake_install_fromsource $version

    # Check version
    if (!(require_new_version $(cmake_version) $version)) {
        cprint "cmake version $(cmake_version) will be used"
    } else {
        if (cmdexists cmake) {
            cerror "cmake version $(cmake_version) will be used but $version is required"
        } else {
            cerror "cmake is not installed"
        }
    }

    cprintend
}