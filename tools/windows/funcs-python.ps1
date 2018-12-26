# 
# Functions related to Python.
# 

. $PSScriptRoot\funcs.ps1
. $PSScriptRoot\funcs-install.ps1
. $PSScriptRoot\funcs-msc.ps1

function python_bin()
{
    if ($global:PYTHONVREQUEST -eq $null) {
        if ((cmdexists python)) {
            return "python"
        } elseif ((cmdexists py)) {
            return "py"
        } else {
            $global:PYTHONVREQUEST=3
        }
    }

    if ((cmdexists py)) {
        return "py -$global:PYTHONVREQUEST"
    } else {
        return "python"
    }
}


function python_full_bin()
{
    return (Get-Command $(python_bin).split()[0]).definition
}


function pip_bin()
{
    if ((python_hasmodule pip)) {
        return "$(python_bin) -m pip"
    } else {
        return "pip"
    }
}


function python_virtualenv_active()
{
    return ((python_get "import sys;print(hasattr(sys, 'real_prefix'))") -eq "True")
}


function python_virtualenv_deactivate()
{
    if ((python_virtualenv_active)) {
        deactivate
    }
}


function python_virtualenv_system_link($arguments)
{
    if ((python_virtualenv_active)) {
        $local:_pkgdir=python_pkg
        $local:_syspkgdir=python_system_pkg
        $local:_restore=Get-Location
        cd $local:_pkgdir
        foreach ($var in $arguments) {
            make-link $local:_syspkgdir/$var $var
        }
        cd $local:_restore
    }
}


function python_hasmodule([string]$module)
{
    $local:restorewd=Get-Location
    cd c:\
    $local:ret=python_get "import sys;sys.stderr=sys.stdout;import $module"
    cd $local:restorewd
    return $local:ret -eq $null
}


function python_get([string]$prog)
{
    if ((python_exists)) {
        Invoke-Expression "$(python_bin) -c ""$prog"""
    }
}


function python_exists()
{
    return cmdexists $(python_bin).split()[0]
}


function python_version()
{
    $local:tmp="import sys;t='{v[0]}.{v[1]}'.format(v=list(sys.version_info[:2]));print(t)"
    $local:tmp=python_get $local:tmp
    if ($local:tmp -eq $null) {
        $local:tmp=-1
    }
    return $local:tmp
}


function python_major_version()
{
    $local:tmp="import sys;print(sys.version_info[0])"
    $local:tmp=python_get $local:tmp 
    if ($local:tmp -eq $null) {
        $local:tmp=-1
    }
    return $local:tmp
}


function python_full_version()
{
    $local:tmp="import sys;t='{v[0]}.{v[1]}.{v[2]}'.format(v=list(sys.version_info[:3]));print(t)"
    $local:tmp=python_get $local:tmp
    if ($local:tmp -eq $null) {
        $local:tmp = -1
    }
    return $local:tmp
}


function python_depdir()
{
    return "dep_$(python_full_version)"
}


function python_include()
{
    return python_get "import distutils.sysconfig; print(distutils.sysconfig.get_python_inc());"
}


function python_lib()
{
    return python_get "import distutils.sysconfig,os; print(os.path.join(distutils.sysconfig.get_config_var('BINDIR'),'python'+distutils.sysconfig.get_config_var('VERSION')+'.dll'));"
}


function python_pkg()
{
    return python_get "import distutils.sysconfig; print(distutils.sysconfig.get_python_lib());"
}


function python_info()
{
    cprint "Python version: $(python_full_version)"
    cprint "Python virtual environment: $(python_virtualenv_active)"
    cprint "Python location: $(python_full_bin)"
    cprint "Python package directory: $(python_pkg)"
    cprint "Python include: $(python_include)"
    cprint "Python library: $(python_lib)"
}


function pip_version()
{
    return (Invoke-Expression "$(pip_bin) --version").split()[1]
}


function pip_source()
{
    return (Invoke-Expression "$(pip_bin) --version").split()[3]
}


function pip_info()
{
    cprint "Pip: $(pip_source)"
}


function pip_install($arguments)
{
    if ((install_systemwide) -or (python_virtualenv_active)) {
        Invoke-Expression "$(pip_bin) install $arguments"
    } else {
        Invoke-Expression "$(pip_bin) install --user $arguments"
    }
}


function pip_upgrade($arguments)
{
    if ((install_systemwide) -or (python_virtualenv_active)) {
        Invoke-Expression "$(pip_bin) install --upgrade $arguments"
    } else {
        Invoke-Expression "$(pip_bin) install --upgrade --user $arguments"
    }
}


function require_python([AllowNull()][string]$version)
{
    if ($version -ne "") {
        $global:PYTHONVREQUEST = $version
    }

    if ($global:PYTHONVREQUEST -eq "") {
        cprint "Verify python (no specific version requested) ..."
    } else {
        cprint "Verify python $global:PYTHONVREQUEST ..."
    }

    # Check version
    if (!(require_new_version_strict $(python_full_version) $global:PYTHONVREQUEST)) {
        cprint "Python version $(python_full_version) is used"
        return
    }

    # Install from source
    python_install_fromsource $global:PYTHONVREQUEST

    # Check version
    if (!(require_new_version_strict $(python_full_version) $global:PYTHONVREQUEST)) {
        cprint "Python version $(python_full_version) is used"
    } else {
        if ((python_exists)) {
            cerror "Python version $(python_full_version) is used but $global:PYTHONVREQUEST is required"
        } else {
            cerror "Python is not installed"
        }
    }
}


function python_install_fromsource([AllowNull()][string]$version)
{
    $local:restorewd = Get-Location

    cprint "Download python ..."
    mkdir python -Force
    cd python

    # Local or latest version that matches the requested version
    $local:lversion = $(get_local_version_strict $version)
    require_web_essentials
    if ($local:lversion -eq $null) {
        $local:version = python_latest $local:version
    } else {
        $local:version = $local:lversion
    }

    # Download the installer
    $local:installerprefix = "python-$local:version"
    $local:installer = Get-ChildItem "$local:installerprefix*.*"
    if ( !(dryrun) -and $local:installer -eq $null ) {
        python_download $local:version
        $local:installer = Get-ChildItem "$local:installerprefix*.*"
    }

    # Run the installer
    if ($local:installer -ne $null) {
        cprint "Install python $local:installer ..."
        if (!(dryrun)) {
            $local:arguments = @{}
            $arguments["msi"] = @()
            $arguments["exe"] = @()

            $arguments["msi"] += "/passive"
            $arguments["exe"] += "/passive"

            $local:systemwide = [int]$(install_systemwide)
            $arguments["msi"] += "ALLUSERS=""$local:systemwide"""
            $arguments["exe"] += "InstallAllUsers=$local:systemwide"

            $arguments["exe"] += "Include_test=0"
            $arguments["exe"] += "PrependPath=1"
            $arguments["exe"] += "Include_launcher=1"
            $arguments["exe"] += "InstallLauncherAllUsers=$local:systemwide"

            install_any $local:installer.Name $local:arguments
            updateBinPath
        }
    } else {
        cerror "Could not find python $global:PYTHONVREQUEST"
    }

    cd $local:restorewd
}


function python_url()
{
    echo "https://www.python.org/ftp/python/"
}


function python_latest($local:rversion) {
    $local:pattern = "(\d\.\d\.\d)"
    if ($local:rversion -ne $null) {
        $local:pattern = $local:rversion.split('.')
        if ($local:pattern.Length -lt 3) {
            $local:pattern += @("\d")*(3-$local:pattern.Length)
        } else {
            $local:pattern = $local:pattern[0..2]
        }
        $local:pattern = [string]::Join("\.",$local:pattern)
        $local:pattern = "($local:pattern)"
    }

    $local:mversion = (0,0,0)
    $local:version = $null
    $local:content = Invoke-WebRequest $(python_url)
    foreach ($link in $local:content.Links) {
        $local:m = [regex]::match($link.href,$local:pattern)
        if ($local:m.Success) {
            if ((python_download_info $m.Groups[1].Value) -eq $null) {
                continue
            }
            $local:tmp = $m.Groups[1].Value.split('.')
            if (($local:tmp.Length) -lt 3) {
                $local:tmp += @(0)*(3-($local:tmp.Length))
            }
            if ($local:tmp[0] -gt $local:mversion[0]) {
                $local:mversion = $local:tmp
                $local:version = $m.Groups[1].Value
            } elseif ($local:tmp[0] -eq $local:mversion[0]) {
                if ($local:tmp[1] -gt $local:mversion[1]) {
                    $local:mversion = $local:tmp
                    $local:version = $m.Groups[1].Value
                } elseif ($local:tmp[1] -eq $local:mversion[1]) {
                    if ($local:tmp[2] -gt $local:mversion[2]) {
                        $local:mversion = $local:tmp
                        $local:version = $m.Groups[1].Value
                    }
                }
            }
        }
    }

    return $local:version
}


function python_download_info([string]$version)
{
    $local:pattern = "(python-$version.[exmsi]{3})"
    if ((install_arch) -eq 64){
        $local:pattern = "(python-$version[-.]amd64.[exmsi]{3})"
    }

    $local:content = Invoke-WebRequest "$(python_url)/$version/"
    foreach ($link in $local:content.Links) {
        $local:m = [regex]::match($link.href,$local:pattern)
        if ($local:m.Success) {
            $local:filename = $m.Groups[1].Value
            return ("$(python_url)/$version/$local:filename","$(Get-Location)\$local:filename")
        }
    }
}

function python_download([string]$version)
{
    $local:arr = python_download_info($version)
    if ($local:arr -ne $null) {
        download_file $local:arr[0] $local:arr[1]
    }
}


function require_pip()
{
    pip_ensure
    pip_upgrade pip
    pip_upgrade setuptools
    pip_upgrade wheel
}


function pip_ensure()
{
    if (!(python_hasmodule pip)) {
        if ((python_hasmodule ensurepip)) {
            if ((python_virtualenv_active)) {
                Invoke-Expression "$(python_bin) -m ensurepip"
            } else {
                if ((install_systemwide)) {
                    Invoke-Expression "$(python_bin) -m ensurepip"
                } else {
                    Invoke-Expression "$(python_bin) -m ensurepip --user"
                }
            }
        }
    }
}


function python_compiler()
{
    $local:s = python_get "import sys;print(sys.version)"
    $local:m = [regex]::match($local:s,"\[(.+?)\]")
    if ($local:m.Success) {
        $local:m = [regex]::match($m.Groups[1].Value,"v\.([\d]+)")
        if ($local:m.Success) {
            return [int]($m.Groups[1].Value)
        }
    }
    
}


function python_arch()
{
    $local:s = python_get "import sys;print(sys.version)"
    $local:m = [regex]::match($local:s,"\[(.+?)\]")
    if ($local:m.Success) {
        $local:m = [regex]::match($m.Groups[1].Value,"64 bit")
        if ($local:m.Success) {
            return 64
        } else {
            return 32
        }
    }
}


function python_init_compiler()
{
    $local:msccompiler = python_compiler
    $local:mscinfo = msc_info $local:msccompiler

    if (!(Test-Path $local:mscinfo["vcvarsall"] -pathType leaf)) {
        cerror "MSC compiler v.$local:msccompiler not installed"
        return
    }

    $local:vcvarsall = $local:mscinfo["vcvarsall"]
    $local:vcvarsall_args = $local:mscinfo[$(python_arch)]
    Invoke-Expression """$local:vcvarsall"" $local:vcvarsall_args"
    cprint $?
}