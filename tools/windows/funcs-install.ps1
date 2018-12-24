#
# Windows installers
#

. $PSScriptRoot\funcs.ps1

# ============install_msi============
# Description: Install a Microsoft installer package
function install_msi($filename,$arguments) {
    
    $msiarguments = @()
    $msiarguments += "/i"
    $msiarguments += "`"$filename`""
    $msiarguments += $arguments
	Start-Process "msiexec.exe" -ArgumentList $msiarguments -Wait
    ThrowIfFailed "Failed to install $filename"
}


# ============install_exe============
# Description: Install an executable
function install_exe($execname,$arguments) {
    if ($arguments.length -gt 0) {
	    Start-Process $execname -ArgumentList $arguments -Wait
    } else {
        Start-Process $execname -Wait
    }
    ThrowIfFailed "Failed to install $execname"
}


# ============install_any============
# Description: Install an package or executable
function install_any($filename,$argdict) {
    if ($filename.endswith(".msi")) {
        install_msi $filename $argdict["msi"]
    } elseif ($filename.endswith(".exe")) {
        install_exe $filename $argdict["exe"]
    } else {
        throw "Don't know how to install $filename"
    }
}

function _addBinPath([string]$path,[string]$target,[bool]$prepend=$false) {
    $envpath = [Environment]::GetEnvironmentVariable("Path", $target)
    if (!($envpath.contains($path))) {
        if ($prepend) {
            [Environment]::SetEnvironmentVariable("Path", "$path;$envpath", $target)
        } else {
            [Environment]::SetEnvironmentVariable("Path", "$envpath;$path", $target)
        }
    }
}

function appendBinPath($path) {
    if ((install_systemwide)) {
        _addBinPath $path "Machine" $false
    } else {
        _addBinPath $path "User" $false
    }
    _addBinPath $path "Process" $false
}

function prependBinPath($path) {
    if ((install_systemwide)) {
        _addBinPath $path "Machine" $true
    } else {
        _addBinPath $path "User" $true
    }
    _addBinPath $path "Process" $true
}

function updateBinPath() {
    $pathmachine = [Environment]::GetEnvironmentVariable("Path", "Machine")
    $pathuser = [Environment]::GetEnvironmentVariable("Path", "User")
    [Environment]::SetEnvironmentVariable("Path", "$pathmachine;$pathuser;$env:path", "Process")
}

function defaultTargetDir() {
    if ((install_systemwide)) {
        if ((install_arch) -eq 64){
            return $env:programfiles
        } else {
            return "$env:programfiles(x86)"
        }
    } else {
        return "$env:localappdata"
    }
}