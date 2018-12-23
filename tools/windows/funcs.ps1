# 
# Helper functions.
# 

. $PSScriptRoot\funcs-string.ps1
. $PSScriptRoot\funcs-essentials.ps1

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
function cprint_error()
{
    Write-Host -ForegroundColor Red $args
}


# ============cprintstart============
# Description: output to stdout in color
# Usage: cprintstart
function cprintstart()
{
    cprint ""
    cprint ""
    cprint ""
    cprint "======================"
}


# ============cprintend============
# Description: output to stdout in color
# Usage: cprintend
function cprintend()
{
    cprint "======================"
    cprint ""
    cprint ""
    cprint ""
}


# ============system_privileges============
# Description: check for root access
function system_privileges()
{  
    $user = [Security.Principal.WindowsIdentity]::GetCurrent()
    (New-Object Security.Principal.WindowsPrincipal $user).IsInRole([Security.Principal.WindowsBuiltinRole]::Administrator)  
}


# ============timer============
# Description: prompt yes/no?
function YesNoQuestion([string]$question)
{
    $local:response = Read-Host "$question ( Y / n ) "
    Switch ($local:response) { 
       Y {return $true} 
       N {return $false} 
       Default {return $true} 
    } 
}


# ============timer============
# Description: start and stop timer
function timer([string]$reset)
{
    $local:seconds = strtoint (Get-Date (Get-Date).ToUniversalTime() -UFormat %s)
    if ($reset -eq "reset") {
        $global:START_TIME = $local:seconds
    } else {
        $local:elapsedtot=$local:seconds - $global:START_TIME
        $local:elapsedmin=[int]($local:elapsedtot/60)
        $local:elapsedsec=[int]($local:elapsedtot-60*$local:elapsedmin)
        cprint "Total execution time = $local:elapsedmin min $local:elapsedsec sec"
    }
}


# ============dryrun============
# Description:
function dryrun([string]$reset,[AllowNull()][boolean]$resetvalue)
{
    if ($global:DRYRUN -eq $null -or $reset -eq "reset") {
        if ($resetvalue -eq $null) {
            $global:DRYRUN = $true
        } else {
            $global:DRYRUN = $resetvalue
        }
        return
    }
    return $global:DRYRUN
}


# ============install_systemwide============
# Description: check for user or system installation
function install_systemwide([string]$reset,[AllowNull()][boolean]$resetvalue)
{
    if ($global:INSTALL_SYSTEMWIDE -eq $null -or $reset -eq "reset") {
        if ($resetvalue -eq $null) {
            $global:INSTALL_SYSTEMWIDE = system_privileges
        } else {
            $global:INSTALL_SYSTEMWIDE = $resetvalue
        }
        if ($reset -eq "reset") {
            return
        }
    }

    if (!(system_privileges)) {
        return $false
    } else {
        return $global:INSTALL_SYSTEMWIDE
    }
}


# ============install_arch============
# Description: installation architecture
function install_arch([string]$reset,[AllowNull()][int]$resetvalue)
{
    if ($global:INSTALL_ARCH -eq $null -or $reset -eq "reset") {
        if ($resetvalue -eq $null) {
            $global:INSTALL_ARCH = -1
        } else {
            $global:INSTALL_ARCH = $resetvalue
        }
        if ($reset -eq "reset") {
            return
        }
    }

    if ($global:INSTALL_ARCH -eq -1) {
        if ([Environment]::Is64BitOperatingSystem) {
            $global:ARG_ARCH = 64
        } else {
            $global:ARG_ARCH = 32
        }
    }

    return $global:ARG_ARCH
}


# ============cmdexists============
# Description: 
function cmdexists([string]$cmd)
{
    return (Get-Command $cmd -errorAction SilentlyContinue) -ne $null
}


# ============project_folder============
# Description: Project folder
function project_folder()
{
    return Resolve-Path $PSScriptRoot\..\..
}


# ============project_name============
# Description: Project name
function project_name()
{
    return (Get-Item (project_folder)).Basename
}


# ============install_info============
# Description: 
function install_info()
{
    cprint "Root priviliges: $(system_privileges)"
    cprint "System wide installation: $(install_systemwide)"
    cprint "Target architecture: $(install_arch)"
    #cprint "Prefix for dependencies: $(project_prefix)"
    #cprint "Opt directory: $(project_opt)"
    #cprint "Resource file: $(project_resource)"
}
