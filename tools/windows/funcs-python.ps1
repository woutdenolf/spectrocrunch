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
    return ((python_get "import sys;print(sys.prefix!=getattr(sys,'base_prefix',sys.prefix))") -eq "True")
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
        $local:_restore = (Get-Location).Path
        cd $local:_pkgdir
        foreach ($var in $arguments) {
            make-link $local:_syspkgdir/$var $var
        }
        cd $local:_restore
    }
}


function python_hasmodule([string]$module)
{
    $local:restorewd = (Get-Location).Path
    cd c:\
    $local:ret=python_get "import sys;sys.stderr=sys.stdout;import $module"
    cd $local:restorewd
    return $local:ret -eq $null
}


function python_exec()
{
    if ((python_exists)) {
        Invoke-Expression "$(python_bin) $args"
    }
}


function python_get([string]$prog)
{
    python_exec -c """$prog"""
}



function python_exists()
{
    return cmdexists $(python_bin).split()[0]
}


function python_minor_version()
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


function python_version()
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
    return "dep_python$(python_version)"
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
    cprint "Python version: $(python_version)"
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
    if (!(require_new_version_strict $(python_version) $global:PYTHONVREQUEST)) {
        cprint "Python version $(python_version) is used"
        return
    }

    # Install from source
    python_install_fromsource $global:PYTHONVREQUEST

    # Check version
    if (!(require_new_version_strict $(python_version) $global:PYTHONVREQUEST)) {
        cprint "Python version $(python_version) is used"
    } else {
        if ((python_exists)) {
            cerror "Python version $(python_version) is used but $global:PYTHONVREQUEST is required"
        } else {
            cerror "Python is not installed"
        }
    }
}


function python_install_fromsource([AllowNull()][string]$version)
{
    $local:restorewd = (Get-Location).Path

    cprint "Download python ..."
    mkdir python -Force
    cd python

    # Local or latest version that matches the requested version
    $local:lversion,$local:installer = get_local_version_strict $version
    if ($local:installer -eq $null) {
        require_web_essentials
        $local:version = python_latest $local:version
        if (!(dryrun))
        {
            $local:installer = python_download $local:version
        }
    } else {
        $local:version = $local:lversion
    }

    # Run the installer
    cprint "Install python $local:version ..."
    if (!(dryrun))
    {
        if ($local:installer -ne $null) {
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
    
                $local:tmp = install_any $local:installer $local:arguments

                # Add to env:path
                updateBinPath
                $local:tmp = version_intarray $local:version 2
                $local:tmp = [string]::Join(".",$local:tmp)
                $local:installpath = registry-value "Registry::HKEY_LOCAL_MACHINE\SOFTWARE\Python\PythonCore\$local:tmp\InstallPath" "(default)"
                if ($local:installpath -ne $null) {
                    prependBinPath "$local:installpath;$local:installpath\Scripts"
                }
            }
        } else {
            cerror "Could not find python $global:PYTHONVREQUEST"
        }
    }
    
    cd $local:restorewd
}


function python_url()
{
    return "https://www.python.org/ftp/python/"
}


function python_latest($local:rversion) {
    $local:pattern = version_pattern $local:rversion 3
    $local:mversion = (0,0,0)
    $local:version = $null
    $local:content = Invoke-WebRequest $(python_url)
    foreach ($link in $local:content.Links) {
        $local:m = [regex]::match($link.href,$local:pattern)
        if ($local:m.Success) {
            if ((python_download_info $m.Groups[1].Value) -eq $null) {
                continue
            }
            $local:tmp = version_intarray $m.Groups[1].Value 3
            if ((version_intarray_higher $local:mversion $local:tmp)) {
                $local:mversion = $local:tmp
                $local:version = $m.Groups[1].Value
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
            return ("$(python_url)/$version/$local:filename",$(joinPath (Get-Location).Path $local:filename))
        }
    }
}

function python_download([string]$version)
{
    $local:arr = python_download_info($version)
    if ($local:arr -ne $null) {
        return download_file $local:arr[0] $local:arr[1]
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


function python_init_compiler([boolean]$strict=$true)
{
    $local:msc_ver = python_compiler
    $local:msc_arch = python_arch
    require_msc $local:msc_ver $local:msc_arch -strict $strict
}
