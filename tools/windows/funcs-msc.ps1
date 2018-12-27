#
# Microsoft Visual C++ compiler
#

. $PSScriptRoot\funcs.ps1

function msc_versions()
{
    $versions = @{}

    $tmp = @{}
    $tmp["version"] = "6.0"
    $tmp["name"] = "Visual Studio 6.0"
    $versions[1200] = $tmp

    $tmp = @{}
    $tmp["version"] = "7.0"
    $tmp["name"] = "Visual Studio .NET 2002"
    $versions[1300] = $tmp

    $tmp = @{}
    $tmp["version"] = "7.1"
    $tmp["name"] = "Visual Studio .NET 2003"
    $versions[1310] = $tmp

    $tmp = @{}
    $tmp["version"] = "8.0"
    $tmp["name"] = "Visual Studio 2005"
    $versions[1400] = $tmp

    # python 2.7
    $tmp = @{}
    $tmp["version"] = "9.0"
    $tmp["name"] = "Visual Studio 2008"
    $tmp["link"] = "https://download.microsoft.com/download/7/9/6/796EF2E4-801B-4FC4-AB28-B59FBF6D907B/VCForPython27.msi"
    $tmp["install_args"] = @()
    $versions[1500] = $tmp

    # python 3.4
    $tmp = @{}
    $tmp["version"] = "10.0"
    $tmp["name"] = "Visual Studio 2010"
    $versions[1600] = $tmp

    $tmp = @{}
    $tmp["version"] = "11.0"
    $tmp["name"] = "Visual Studio 2012"
    $versions[1700] = $tmp

    $tmp = @{}
    $tmp["version"] = "12.0"
    $tmp["name"] = "Visual Studio 2013"
    $versions[1800] = $tmp

    # python 3.5, 3.6
    $tmp = @{}
    $tmp["version"] = "14.0"
    $tmp["name"] = "Visual Studio 2015"
    $tmp["vcvarsall"] = "C:\Program Files (x86)\Microsoft Visual Studio\2017\BuildTools\VC\Auxiliary\Build\vcvarsall.bat"
    $tmp[32] = "x86 -vcvars_ver=14.0"
    $tmp[64] = "x64 -vcvars_ver=14.0"
    $tmp["link"] = "https://www.visualstudio.com/thank-you-downloading-visual-studio/?sku=BuildTools&rel=15#"
    $tmp["install_args"] = @()
    $tmp["arguments"] += "--add Microsoft.VisualStudio.Component.VC.140"
    $versions[1900] = $tmp

    $tmp = @{}
    $tmp["version"] = "15.0"
    $tmp["name"] = "Visual Studio 2017"
    $tmp["vcvarsall"] = "C:\Program Files (x86)\Microsoft Visual Studio\2017\BuildTools\VC\Auxiliary\Build\vcvarsall.bat"
    $tmp[32] = "x86 -vcvars_ver=15.0"
    $tmp[64] = "x64 -vcvars_ver=15.0"
    $tmp["link"] = "https://www.visualstudio.com/thank-you-downloading-visual-studio/?sku=BuildTools&rel=15#"
    $tmp["install_args"] = @()
    $versions[1910] = $tmp

    $tmp = @{}
    $tmp["version"] = "15.3"
    $tmp["name"] = "Visual Studio 2017"
    $tmp["vcvarsall"] = "C:\Program Files (x86)\Microsoft Visual Studio\2017\BuildTools\VC\Auxiliary\Build\vcvarsall.bat"
    $tmp[32] = "x86 -vcvars_ver=15.3"
    $tmp[64] = "x64 -vcvars_ver=15.3"
    $tmp["link"] = "https://www.visualstudio.com/thank-you-downloading-visual-studio/?sku=BuildTools&rel=15#"
    $tmp["install_args"] = @()
    $versions[1911] = $tmp

    $tmp = @{}
    $tmp["version"] = "15.5"
    $tmp["name"] = "Visual Studio 2017"
    $tmp["vcvarsall"] = "C:\Program Files (x86)\Microsoft Visual Studio\2017\BuildTools\VC\Auxiliary\Build\vcvarsall.bat"
    $tmp[32] = "x86 -vcvars_ver=15.5"
    $tmp[64] = "x64 -vcvars_ver=15.5"
    $tmp["link"] = "https://www.visualstudio.com/thank-you-downloading-visual-studio/?sku=BuildTools&rel=15#"
    $tmp["install_args"] = @()
    $versions[1912] = $tmp

    $tmp = @{}
    $tmp["version"] = "15.6"
    $tmp["name"] = "Visual Studio 2017"
    $tmp["vcvarsall"] = "C:\Program Files (x86)\Microsoft Visual Studio\2017\BuildTools\VC\Auxiliary\Build\vcvarsall.bat"
    $tmp[32] = "x86 -vcvars_ver=15.6"
    $tmp[64] = "x64 -vcvars_ver=15.6"
    $tmp["link"] = "https://www.visualstudio.com/thank-you-downloading-visual-studio/?sku=BuildTools&rel=15#"
    $tmp["install_args"] = @()
    $versions[1913] = $tmp

    $tmp = @{}
    $tmp["version"] = "15.7"
    $tmp["name"] = "Visual Studio 2017"
    $tmp["vcvarsall"] = "C:\Program Files (x86)\Microsoft Visual Studio\2017\BuildTools\VC\Auxiliary\Build\vcvarsall.bat"
    $tmp[32] = "x86 -vcvars_ver=15.7"
    $tmp[64] = "x64 -vcvars_ver=15.7"
    $tmp["link"] = "https://www.visualstudio.com/thank-you-downloading-visual-studio/?sku=BuildTools&rel=15#"
    $tmp["install_args"] = @()
    $versions[1914] = $tmp

    $tmp = @{}
    $tmp["version"] = "15.8"
    $tmp["name"] = "Visual Studio 2017"
    $tmp["vcvarsall"] = "C:\Program Files (x86)\Microsoft Visual Studio\2017\BuildTools\VC\Auxiliary\Build\vcvarsall.bat"
    $tmp[32] = "x86 -vcvars_ver=15.8"
    $tmp[64] = "x64 -vcvars_ver=15.8"
    $tmp["link"] = "https://www.visualstudio.com/thank-you-downloading-visual-studio/?sku=BuildTools&rel=15#"
    $tmp["install_args"] = @()
    $versions[1915] = $tmp

    # python 3.7
    $tmp = @{}
    $tmp["version"] = "15.9"
    $tmp["name"] = "Visual Studio 2017"
    $tmp["vcvarsall"] = "C:\Program Files (x86)\Microsoft Visual Studio\2017\BuildTools\VC\Auxiliary\Build\vcvarsall.bat"
    $tmp[32] = "x86 -vcvars_ver=15.9"
    $tmp[64] = "x64 -vcvars_ver=15.9"
    $tmp["link"] = "https://www.visualstudio.com/thank-you-downloading-visual-studio/?sku=BuildTools&rel=15#"
    $tmp["install_args"] = @()
    $versions[1916] = $tmp

    return $versions
}

function msc_info([int]$msc_ver)
{
    return $(msc_versions)[$msc_ver]
}

function require_msc([int]$msc_ver,[int]$arch)
{
    $local:mscinfo = msc_info $msc_ver
    if (Test-Path $local:mscinfo["vcvarsall"] -pathType leaf) {
        install_msc $local:mscinfo
    }

    if (Test-Path $local:mscinfo["vcvarsall"] -pathType leaf) {
        cerror "$local:mscinfo[""name""] ($local:mscinfo[""version""], $msc_ver) not installed"
    } else {
        $local:vcvarsall = $local:mscinfo["vcvarsall"]
        $local:vcvarsall_args = $local:mscinfo[$arch]
        Invoke-Expression """$local:vcvarsall"" $local:vcvarsall_args"
        if ($?) {
            cerror "$local:mscinfo[""name""] ($local:mscinfo[""version""], $msc_ver) not installed"
        }
    }
}

function install_msc($mscinfo)
{
    vs_buildtools.exe --update --quiet --wait
    vs_enterprise.exe update --wait --passive --norestart
}