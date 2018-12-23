# 
# Functions related to Python.
# 

. $PSScriptRoot\funcs.ps1

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
    if ((cmdexists $(python_bin).split()[0])) {
        Invoke-Expression "$(python_bin) -c ""$prog"""
    }
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
        $local:tmp=-1
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
    python_get "import distutils.sysconfig; print(distutils.sysconfig.get_python_lib());"
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