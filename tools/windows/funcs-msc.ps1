#
# Microsoft Visual C++ compiler
#

. $PSScriptRoot\funcs.ps1
. $PSScriptRoot\funcs-install.ps1

function require_vssetup()
{
    cprint "Checking VSSetup ..."

    if (cmdexists Get-VSSetupInstance) {
        cprint "VSSetup is installed"
        return
    }

    cprint "Install VSSetup ..."
    if (!(dryrun)) {
        require_web_essentials

        if (cmdexists Install-Module) {
            Install-Module VSSetup -Scope CurrentUser
        } else {
            $local:filename = download_git_release "Microsoft" "vssetup.powershell" ".zip"
            if ($local:filename -ne $null) {
                Unzip $local:filename "$([Environment]::GetFolderPath("MyDocuments"))\WindowsPowerShell\Modules\VSSetup"
            }
        }
    
        if (cmdexists Get-VSSetupInstance) {
            cprint "VSSetup is installed"
        } else {
            cerror "VSSetup is not installed"
        }
    }
}

function msc_versions()
{
    $versions = @{}

    $tmp = @{}
    $tmp["msc_ver"] = 1200
    $tmp["version"] = "6.0"
    $tmp["name"] = "Visual Studio 6.0"
    $versions[$tmp["msc_ver"]] = $tmp

    $tmp = @{}
    $tmp["msc_ver"] = 1300
    $tmp["version"] = "7.0"
    $tmp["name"] = "Visual Studio .NET 2002"
    $versions[$tmp["msc_ver"]] = $tmp

    $tmp = @{}
    $tmp["msc_ver"] = 1310
    $tmp["version"] = "7.1"
    $tmp["name"] = "Visual Studio .NET 2003"
    $versions[$tmp["msc_ver"]] = $tmp

    $tmp = @{}
    $tmp["msc_ver"] = 1400
    $tmp["version"] = "8.0"
    $tmp["name"] = "Visual Studio 2005"
    $versions[$tmp["msc_ver"]] = $tmp

    # python 2.7
    $tmp = @{}
    $tmp["msc_ver"] = 1500
    $tmp["version"] = "9.0"
    $tmp["name"] = "Visual Studio 2008"
    $tmp["cmake"] = @{}
    $tmp["cmake"][32] = @{}
    $tmp["cmake"][64] = @{}
    $tmp["cmake"][32]["generator"] = "Visual Studio 9 2008"
    $tmp["cmake"][32]["arch"] = "Win32"
    $tmp["cmake"][32]["toolset"] = $null
    $tmp["cmake"][64]["generator"] = "Visual Studio 9 2008"
    $tmp["cmake"][64]["arch"] = "x64"
    $tmp["cmake"][64]["toolset"] = $null
    $tmp["vcvarsall"] = "Programs\Common\Microsoft\Visual C++ for Python\9.0\vcvarsall.bat"
    $tmp[32] = "x86"
    $tmp[64] = "amd64"
    $tmp["link"] = "https://download.microsoft.com/download/7/9/6/796EF2E4-801B-4FC4-AB28-B59FBF6D907B/VCForPython27.msi"
    $tmp["filename"] = "VCForPython27.msi"
    $tmp["install_args"] = @()
    $tmp["install_args"] += "/quiet"
    $tmp["install_args"] += "/passive"
    $local:systemwide = [int]$(install_systemwide)
    $tmp["install_args"] += "ALLUSERS=""$local:systemwide"""
    $versions[$tmp["msc_ver"]] = $tmp

    # python 3.4
    $tmp = @{}
    $tmp["msc_ver"] = 1600
    $tmp["version"] = "10.0"
    $tmp["name"] = "Visual Studio 2010"
    $tmp["cmake"] = @{}
    $tmp["cmake"][32] = @{}
    $tmp["cmake"][64] = @{}
    $tmp["cmake"][32]["generator"] = "Visual Studio 10 2010"
    $tmp["cmake"][32]["arch"] = "Win32"
    $tmp["cmake"][32]["toolset"] = $null
    $tmp["cmake"][64]["generator"] = "Visual Studio 10 2010"
    $tmp["cmake"][64]["arch"] = "x64"
    $tmp["cmake"][64]["toolset"] = $null
    $tmp["vcvarsall"] = "Microsoft Visual Studio 10.0\VC\vcvarsall.bat"
    $tmp[32] = "x86"
    $tmp[64] = "amd64"
    $versions[$tmp["msc_ver"]] = $tmp

    $tmp = @{}
    $tmp["msc_ver"] = 1700
    $tmp["version"] = "11.0"
    $tmp["name"] = "Visual Studio 2012"
    $tmp["cmake"] = @{}
    $tmp["cmake"][32] = @{}
    $tmp["cmake"][64] = @{}
    $tmp["cmake"][32]["generator"] = "Visual Studio 11 2012"
    $tmp["cmake"][32]["arch"] = "Win32"
    $tmp["cmake"][32]["toolset"] = $null
    $tmp["cmake"][64]["generator"] = "Visual Studio 11 2012"
    $tmp["cmake"][64]["arch"] = "x64"
    $tmp["cmake"][64]["toolset"] = $null
    $tmp["vcvarsall"] = "Microsoft Visual Studio 11.0\VC\vcvarsall.bat"
    $versions[$tmp["msc_ver"]] = $tmp

    $tmp = @{}
    $tmp["msc_ver"] = 1800
    $tmp["version"] = "12.0"
    $tmp["name"] = "Visual Studio 2013"
    $tmp["cmake"] = @{}
    $tmp["cmake"][32] = @{}
    $tmp["cmake"][64] = @{}
    $tmp["cmake"][32]["generator"] = "Visual Studio 12 2013"
    $tmp["cmake"][32]["arch"] = "Win32"
    $tmp["cmake"][32]["toolset"] = $null
    $tmp["cmake"][64]["generator"] = "Visual Studio 12 2013"
    $tmp["cmake"][64]["arch"] = "x64"
    $tmp["cmake"][64]["toolset"] = $null
    $tmp["vcvarsall"] = "Microsoft Visual Studio 12.0\VC\vcvarsall.bat"
    $versions[$tmp["msc_ver"]] = $tmp

    #https://docs.microsoft.com/en-us/visualstudio/install/workload-component-id-vs-build-tools?view=vs-2017
    $vs_buildtools = "https://download.visualstudio.microsoft.com/download/pr/a46d2db7-bd7b-43ee-bd7b-12624297e4ec/11b9c9bd44ec2b475f6da3d1802b3d00/vs_buildtools.exe"

    # python 3.5
    $tmp = @{}
    $tmp["msc_ver"] = 1900
    $tmp["version"] = "14.0"
    $tmp["name"] = "Visual Studio 2015"
    $tmp["cmake"] = @{}
    $tmp["cmake"][32] = @{}
    $tmp["cmake"][64] = @{}
    $tmp["cmake"][32]["generator"] = "Visual Studio 14 2015"
    $tmp["cmake"][32]["arch"] = "Win32"
    $tmp["cmake"][32]["toolset"] = "140"
    $tmp["cmake"][64]["generator"] = "Visual Studio 14 2015"
    $tmp["cmake"][64]["arch"] = "x64"
    $tmp["cmake"][64]["toolset"] = "140,host=x64"
    $tmp["manager"] = "vssetup"
    $tmp["installer"] = "Microsoft Visual Studio\Installer\vs_installershell.exe"
    $tmp["toolset"] = "Microsoft.VisualStudio.Component.VC.140"
    $tmp["vcvarsall"] = "Microsoft Visual Studio 14.0\VC\vcvarsall.bat"
    $tmp[32] = "x86"
    $tmp[64] = "x64"
    $tmp["link"] = $vs_buildtools
    $tmp["filename"] = "vs_buildtools.exe"
    $tmp["install_args"] = @()
    $tmp["install_args"] += "--quiet"
    $versions[$tmp["msc_ver"]] = $tmp

    $tmp = @{}
    $tmp["msc_ver"] = 1910
    $tmp["version"] = "15.0"
    $tmp["name"] = "Visual Studio 2017"
    $tmp["cmake"] = @{}
    $tmp["cmake"][32] = @{}
    $tmp["cmake"][64] = @{}
    $tmp["cmake"][32]["generator"] = "Visual Studio 15 2017"
    $tmp["cmake"][32]["arch"] = "Win32"
    $tmp["cmake"][32]["toolset"] = $null
    $tmp["cmake"][64]["generator"] = "Visual Studio 15 2017"
    $tmp["cmake"][64]["arch"] = "x64"
    $tmp["cmake"][64]["toolset"] = $null
    $versions[$tmp["msc_ver"]] = $tmp

    $tmp = @{}
    $tmp["msc_ver"] = 1911
    $tmp["version"] = "15.3"
    $tmp["name"] = "Visual Studio 2017"
    $tmp["cmake"] = @{}
    $tmp["cmake"][32] = @{}
    $tmp["cmake"][64] = @{}
    $tmp["cmake"][32]["generator"] = "Visual Studio 15 2017"
    $tmp["cmake"][32]["arch"] = "Win32"
    $tmp["cmake"][32]["toolset"] = "version=14.11"
    $tmp["cmake"][64]["generator"] = "Visual Studio 15 2017"
    $tmp["cmake"][64]["arch"] = "x64"
    $tmp["cmake"][64]["toolset"] = "host=x64,version=14.11"
    $tmp["manager"] = "vssetup"
    $tmp["installer"] = "Microsoft Visual Studio\Installer\vs_installershell.exe"
    $tmp["toolset"] = "Microsoft.VisualStudio.Component.VC.Tools.14.11"
    $tmp["vcvarsall"] = "VC\Auxiliary\Build\vcvarsall.bat"
    $tmp[32] = "x86 -vcvars_ver=14.11"
    $tmp[64] = "x64 -vcvars_ver=14.11"
    $tmp["link"] = $vs_buildtools
    $tmp["filename"] = "vs_buildtools.exe"
    $tmp["install_args"] = @()
    $tmp["install_args"] += "--quiet"
    $versions[$tmp["msc_ver"]] = $tmp

    $tmp = @{}
    $tmp["msc_ver"] = 1912
    $tmp["version"] = "15.5"
    $tmp["name"] = "Visual Studio 2017"
    $tmp["cmake"] = @{}
    $tmp["cmake"][32] = @{}
    $tmp["cmake"][64] = @{}
    $tmp["cmake"][32]["generator"] = "Visual Studio 15 2017"
    $tmp["cmake"][32]["arch"] = "Win32"
    $tmp["cmake"][32]["toolset"] = "version=14.12"
    $tmp["cmake"][64]["generator"] = "Visual Studio 15 2017"
    $tmp["cmake"][64]["arch"] = "x64"
    $tmp["cmake"][64]["toolset"] = "host=x64,version=14.12"
    $tmp["manager"] = "vssetup"
    $tmp["installer"] = "Microsoft Visual Studio\Installer\vs_installershell.exe"
    $tmp["toolset"] = "Microsoft.VisualStudio.Component.VC.Tools.14.12"
    $tmp["vcvarsall"] = "VC\Auxiliary\Build\vcvarsall.bat"
    $tmp[32] = "x86 -vcvars_ver=14.12"
    $tmp[64] = "x64 -vcvars_ver=14.12"
    $tmp["link"] = $vs_buildtools
    $tmp["filename"] = "vs_buildtools.exe"
    $tmp["install_args"] = @()
    $tmp["install_args"] += "--quiet"
    $versions[$tmp["msc_ver"]] = $tmp

    $tmp = @{}
    $tmp["msc_ver"] = 1913
    $tmp["version"] = "15.6"
    $tmp["name"] = "Visual Studio 2017"
    $tmp["cmake"] = @{}
    $tmp["cmake"][32] = @{}
    $tmp["cmake"][64] = @{}
    $tmp["cmake"][32]["generator"] = "Visual Studio 15 2017"
    $tmp["cmake"][32]["arch"] = "Win32"
    $tmp["cmake"][32]["toolset"] = "version=14.13"
    $tmp["cmake"][64]["generator"] = "Visual Studio 15 2017"
    $tmp["cmake"][64]["arch"] = "x64"
    $tmp["cmake"][64]["toolset"] = "host=x64,version=14.13"
    $tmp["manager"] = "vssetup"
    $tmp["installer"] = "Microsoft Visual Studio\Installer\vs_installershell.exe"
    $tmp["toolset"] = "Microsoft.VisualStudio.Component.VC.Tools.14.13"
    $tmp["vcvarsall"] = "VC\Auxiliary\Build\vcvarsall.bat"
    $tmp[32] = "x86 -vcvars_ver=14.13"
    $tmp[64] = "x64 -vcvars_ver=14.13"
    $tmp["link"] = $vs_buildtools
    $tmp["filename"] = "vs_buildtools.exe"
    $tmp["install_args"] = @()
    $tmp["install_args"] += "--quiet"
    $versions[$tmp["msc_ver"]] = $tmp

    $tmp = @{}
    $tmp["msc_ver"] = 1914
    $tmp["version"] = "15.7"
    $tmp["name"] = "Visual Studio 2017"
    $tmp["cmake"] = @{}
    $tmp["cmake"][32] = @{}
    $tmp["cmake"][64] = @{}
    $tmp["cmake"][32]["generator"] = "Visual Studio 15 2017"
    $tmp["cmake"][32]["arch"] = "Win32"
    $tmp["cmake"][32]["toolset"] = "version=14.14"
    $tmp["cmake"][64]["generator"] = "Visual Studio 15 2017"
    $tmp["cmake"][64]["arch"] = "x64"
    $tmp["cmake"][64]["toolset"] = "host=x64,version=14.14"
    $tmp["manager"] = "vssetup"
    $tmp["installer"] = "Microsoft Visual Studio\Installer\vs_installershell.exe"
    $tmp["toolset"] = "Microsoft.VisualStudio.Component.VC.Tools.14.14"
    $tmp["vcvarsall"] = "VC\Auxiliary\Build\vcvarsall.bat"
    $tmp[32] = "x86 -vcvars_ver=14.14"
    $tmp[64] = "x64 -vcvars_ver=14.14"
    $tmp["link"] = $vs_buildtools
    $tmp["filename"] = "vs_buildtools.exe"
    $tmp["install_args"] = @()
    $tmp["install_args"] += "--quiet"
    $versions[$tmp["msc_ver"]] = $tmp

    $tmp = @{}
    $tmp["msc_ver"] = 1915
    $tmp["version"] = "15.8"
    $tmp["name"] = "Visual Studio 2017"
    $tmp["cmake"] = @{}
    $tmp["cmake"][32] = @{}
    $tmp["cmake"][64] = @{}
    $tmp["cmake"][32]["generator"] = "Visual Studio 15 2017"
    $tmp["cmake"][32]["arch"] = "Win32"
    $tmp["cmake"][32]["toolset"] = "version=14.15"
    $tmp["cmake"][64]["generator"] = "Visual Studio 15 2017"
    $tmp["cmake"][64]["arch"] = "x64"
    $tmp["cmake"][64]["toolset"] = "host=x64,version=14.15"
    $tmp["manager"] = "vssetup"
    $tmp["installer"] = "Microsoft Visual Studio\Installer\vs_installershell.exe"
    $tmp["toolset"] = "Microsoft.VisualStudio.Component.VC.Tools.14.15"
    $tmp["vcvarsall"] = "VC\Auxiliary\Build\vcvarsall.bat"
    $tmp[32] = "x86 -vcvars_ver=14.15"
    $tmp[64] = "x64 -vcvars_ver=14.15"
    $tmp["link"] = $vs_buildtools
    $tmp["filename"] = "vs_buildtools.exe"
    $tmp["install_args"] = @()
    $tmp["install_args"] += "--quiet"
    $versions[$tmp["msc_ver"]] = $tmp

    # python 3.6,3.7
    $tmp = @{}
    $tmp["msc_ver"] = 1916
    $tmp["version"] = "15.9"
    $tmp["name"] = "Visual Studio 2017"
    $tmp["cmake"] = @{}
    $tmp["cmake"][32] = @{}
    $tmp["cmake"][64] = @{}
    $tmp["cmake"][32]["generator"] = "Visual Studio 15 2017"
    $tmp["cmake"][32]["arch"] = "Win32"
    $tmp["cmake"][32]["toolset"] = $null
    $tmp["cmake"][64]["generator"] = "Visual Studio 15 2017"
    $tmp["cmake"][64]["arch"] = "x64"
    $tmp["cmake"][64]["toolset"] = "host=x64"
    $tmp["manager"] = "vssetup"
    $tmp["installer"] = "Microsoft Visual Studio\Installer\vs_installershell.exe"
    $tmp["toolset"] = "Microsoft.VisualStudio.Component.VC.Tools.x86.x64"
    $tmp["vcvarsall"] = "VC\Auxiliary\Build\vcvarsall.bat"
    $tmp[32] = "x86"
    $tmp[64] = "x64"
    $tmp["link"] = $vs_buildtools
    $tmp["filename"] = "vs_buildtools.exe"
    $tmp["install_args"] = @()
    $tmp["install_args"] += "--quiet"
    $versions[$tmp["msc_ver"]] = $tmp

    return $versions
}

function msc_info([int]$msc_ver)
{
    return $(msc_versions)[$msc_ver]
}


function require_msc([int]$msc_ver,[int]$msc_arch,[boolean]$strict=$true)
{
    if (!$strict) {
        $msc_ver_ = current_msc_ver
        if ($msc_ver_ -ne 0) {
            cprint "MSC version $msc_ver_ is already initialized"
            return
        }
    }

    # Check whether msc_ver is known
    $local:mscinfo = msc_info $msc_ver
    if ($local:mscinfo -eq $null) {
        cerror "MSC version $msc_ver is unknown"
    }

    # Compiler is already in environment
    if ((current_msc_ver) -eq $msc_ver) {
        cprint "MSC version $msc_ver is already initialized"
        return
    }

    # Make sure the compiler is installed
    ensure_msc $msc_ver $msc_arch

    # Initialize the compiler env
    init_msc $msc_ver $msc_arch -force=$false
}


function init_msc([int]$msc_ver,[int]$msc_arch,[boolean]$force=$true)
{
    # Check whether msc_ver is known
    $local:mscinfo = msc_info $msc_ver
    if ($local:mscinfo -eq $null) {
        cerror "MSC version $msc_ver is unknown"
        return
    }

    # Compiler is already in environment
    if ((current_msc_ver) -eq $msc_ver -and !$force) {
        cprint "MSC version $msc_ver is initialized"
        return
    }

    # Initialize environment
    $local:vcvarsall = get_vcvarsall $local:mscinfo
    if ($local:vcvarsall -ne $null) {
        cprint "Initializing MSC environment..."
        # TODO: this does not add the path to cl.exe to the 
        $local:errorcode = Invoke-CmdScript "$local:vcvarsall" $local:mscinfo[$msc_arch]
        if ($local:errorcode -ne 0) {
            cerror "MSC initialization failed"
        }
    }

    # Check compiler version
    if ((current_msc_ver) -eq $msc_ver) {
        cprint "MSC $msc_ver is initialized"
    } else {
        cerror "MSC $msc_ver is not initialized"
    }
}


function ensure_msc([int]$msc_ver,[int]$msc_arch)
{
    # Check whether we know about this version
    $local:mscinfo = msc_info $msc_ver
    if ($local:mscinfo -eq $null) {
        cerror "MSC version $msc_ver is unknown"
        return
    }

    # Install
    if ($local:mscinfo["manager"] -eq "vssetup") {
        require_msc_vssetup $local:mscinfo
    } else {
        require_msc_default $local:mscinfo
    }
}


function current_msc_ver()
{
    $local:clexe = @()
    if (cmdexists cl) {
        $local:clexe += "cl"
    } else {
        if (cmdexists msbuild) {
            $local:vsmanager = Get-VSSetupInstance
            if ($local:vsmanager -ne $null) {
                $local:include = joinPath $local:vsmanager.InstallationPath "VC\Tools\MSVC\*\bin\*\*\cl.exe"
                $local:clexe += Get-ChildItem -Path $local:include | % {"""$_"""}
            }
        }
    }

    foreach ($name in $local:clexe) {
        $local:tmp = & cmd /c "$name 2>&1"
        $local:m = [regex]::match($local:tmp,"([\.\d]+)")
        if ($local:m.Success) {
            return [int][string]::Join("",$m.Groups[1].Value.split('.')[0..1])
        }
    }

    return 0
}


function current_msc_arch()
{
    $local:clexe = @()
    if (cmdexists cl) {
        $local:clexe += "cl"
    } else {
        if (cmdexists msbuild) {
            $local:vsmanager = Get-VSSetupInstance
            if ($local:vsmanager -ne $null) {
                $local:include = joinPath $local:vsmanager.InstallationPath "VC\Tools\MSVC\*\bin\*\*\cl.exe"
                $local:clexe += Get-ChildItem -Path $local:include | % {"""$_"""}
            }
        }
    }

    foreach ($name in $local:clexe) {
        $local:tmp = & cmd /c "$name 2>&1"
        $local:m = [regex]::match($local:tmp,"for x64")
        if ($local:m.Success) {
            return 64
        } else {
            return 32
        }
    }

    return install_arch
}

function require_msc_vssetup($mscinfo)
{
    # Make sure the Visual Studio installer is installed
    require_vsinstaller $local:mscinfo

    # Check Visual Studio installer
    $local:vsinstaller = Get-VSSetupInstance
    if ($local:vsinstaller -eq $null) {
        cerror "Visual Studio installer not installed"
        return
    } else {
        vscomponents_list
    }
    
    # Make sure the requested compiler toolset is installed
    require_vstoolset $local:mscinfo
}


function require_msc_default($mscinfo)
{
    $local:vcvarsall = get_vcvarsall $local:mscinfo
    if ($local:vcvarsall -ne $null) {
        cprint "MSC is already installed..."
        return
    }

    if ($local:mscinfo["link"] -eq $null) {
        cerror "No download link specified for MSC version $($local:mscinfo["$msc_ver"])"
        return
    }
    
    cprint "Download and install MSC ..."
    if (!(dryrun)) {
        $local:filename = joinPath (Get-Location).Path $local:mscinfo["filename"]
        $local:filename = download_file $local:mscinfo["link"] $local:filename
        install_msi $local:filename $local:mscinfo["install_args"]
    }
}


function require_vsinstaller($mscinfo)
{
    require_vssetup

    if ((Get-VSSetupInstance) -eq $null) {
        $local:filename = joinPath (Get-Location).Path $mscinfo["filename"]
        
        cprint "Download Visual Studio installer ..."
        if (!(dryrun)) {
            require_web_essentials
            $local:filename = download_file $mscinfo["link"] $local:filename
        }

        cprint "Install Visual Studio installer ..."
        if (!(dryrun)) {
            cprint $mscinfo["install_args"]
            $local:tmp = install_exe $local:filename $mscinfo["install_args"]
        }
    }
}


function require_vstoolset($mscinfo)
{
    require_component $mscinfo["toolset"]
}


function require_sdk($mscinfo)
{
    require_component $mscinfo["sdk"]
}


function require_component($component)
{
    if (!(vscomponent_installed $local:component)) {
        $local:vsmanager = Get-VSSetupInstance

        cprint "Install component $local:component ..."
        $local:installer = vsinstaller $mscinfo
        if ($local:installer -eq $null) {
            cerror "Visual studio installer not found"
        } else {
            if (!(dryrun)) {
                require_web_essentials
                $local:args = @()
                $local:args += "modify"
                $local:args += "--installPath ""$($vsmanager.InstallationPath)"""
                $local:args += "--add $local:component"
                #$local:args += "--includeRecommended"
                $local:args += "--passive"
                $local:tmp = install_exe $local:installer $local:args $true
            }
        }
    }
    
    if ((vscomponent_installed $local:component)) {
        cprint "MSC component $local:component is installed"
    } else {
        cerror "MSC component $local:component is not installed"
    }
}


function vsinstaller($mscinfo)
{
    $local:bases = ($env:programfiles,${env:programfiles(x86)},$env:localappdata)
    foreach ($base in $local:bases) {
        $local:filename = joinPath $base $mscinfo["installer"]
        if ((Test-Path $local:filename -pathType leaf)) {
            return $local:filename
        }
    }
}


function vscomponents_list()
{
    foreach ($pkg in (Get-VSSetupInstance).packages) {
        if ($pkg.type -in 'workload', 'component') {
            cprint $pkg.type ":  " $pkg.Id
        }
    }

}


function vscomponent_installed([string]$component)
{
    foreach ($pkg in (Get-VSSetupInstance).packages) {
        if ($pkg.Id -eq $component) {
            return $true
        }
    }
    return $false
}


function get_vcvarsall($mscinfo)
{
    if ($local:mscinfo["manager"] -eq "vssetup") {
        $local:vsmanager = Get-VSSetupInstance
        if ($local:vsmanager -ne $null) {
            $local:filename = joinPath $local:vsmanager.InstallationPath $local:mscinfo["vcvarsall"]
            if ((file_exists $local:filename)) {
                cprint "vcvarsall: $local:filename"
                return $local:filename
            }
        }
    }
    
    if ($local:mscinfo["vcvarsall"] -ne $null) {
        $local:bases = ($env:programfiles,${env:programfiles(x86)},$env:localappdata)
        foreach ($base in $local:bases) {
            $local:filename = joinPath $base $local:mscinfo["vcvarsall"]
            if (file_exists $local:filename) {
                cprint "vcvarsall: $local:filename"
                return $local:filename
            }
        }
    }

    cerror "vcvarsall: " $local:mscinfo["vcvarsall"] "Not Found"
}


function msc_cmake_info()
{
    $local:info = @{}
    $local:msc_ver = current_msc_ver
    $local:mscinfo = msc_info $local:msc_ver
    return $mscinfo["cmake"][$(current_msc_arch)]
}