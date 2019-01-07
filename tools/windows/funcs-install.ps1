#
# Windows installers
#

. $PSScriptRoot\funcs.ps1

# ============run_new_process============
# Description:
function run_in_new_process($execname,$arguments,[AllowNull()][boolean]$asadmin) {
    if ($asadmin -eq $true) {
        if ($arguments.length -gt 0) {
            $process = Start-Process "$execname" -ArgumentList $arguments -PassThru -Verb runas
        } else {
            $process = Start-Process "$execname" -PassThru -Verb runas
        }
        while ($process.HandleCount -ne $null) {
            $tmp = Get-Process -InputObject $process
            sleep 1
        }
    } else {
        if ($arguments.length -gt 0) {
            $process = Start-Process "$execname" -ArgumentList $arguments -PassThru
        } else {
            $process = Start-Process "$execname" -PassThru
        }
        Wait-Process -InputObject $process
    }
    if ($process.ExitCode -ne 0 -and $process.ExitCode -ne $null) {
        cerror "Process running ""$execname"" failed"
    }
    return $process.ExitCode
}


# ============install_msi============
# Description: Install a Microsoft installer package
function install_msi($filename,$arguments,[AllowNull()][boolean]$asadmin) {
    
    $msiarguments = @()
    $msiarguments += "/i"
    $msiarguments += """$filename"""
    $msiarguments += $arguments
    run_in_new_process "msiexec.exe" $msiarguments $asadmin
}


# ============install_exe============
# Description: Install an executable
function install_exe($filename,$arguments,[AllowNull()][boolean]$asadmin) {
    run_in_new_process $filename $arguments $asadmin
}


# ============install_any============
# Description: Install an package or executable
function install_any($filename,$argdict,[AllowNull()][boolean]$asadmin) {
    if ($filename.endswith(".msi")) {
        install_msi $filename $argdict["msi"] $asadmin
    } elseif ($filename.endswith(".exe")) {
        install_exe $filename $argdict["exe"] $asadmin
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


function _cleanBinPath([string]$target) {
    $envpath = [Environment]::GetEnvironmentVariable("Path", $target).split(";")
    $newpath = @()
    foreach ($p in $envpath) {
        $pfull = [System.Environment]::ExpandEnvironmentVariables($p)
        if ((Test-Path -Path $pfull)){
            $newpath += $p
        } else {
            cprint "Removing: $p"
        }
    }
    [Environment]::SetEnvironmentVariable("Path", [string]::join(";",$newpath), $target)
}



# ============appendBinPath============
# Description: add to user/machine + process environment
function appendBinPath($path) {
    if ((install_systemwide)) {
        _addBinPath $path "Machine" $false
    } else {
        _addBinPath $path "User" $false
    }
    _addBinPath $path "Process" $false
}


# ============prependBinPath============
# Description: add to user/machine + process environment
function prependBinPath($path) {
    if ((install_systemwide)) {
        _addBinPath $path "Machine" $true
    } else {
        _addBinPath $path "User" $true
    }
    _addBinPath $path "Process" $true
}


# ============cleanBinPath============
# Description: remove non-existing paths from user/machine + process environment
function cleanBinPath() {
    if ((install_systemwide)) {
        _cleanBinPath "Machine"
    } else {
        _cleanBinPath "User"
    }
    _cleanBinPath "Process"
}


# ============updateBinPath============
# Description: update path environement variable
function updateBinPath() {
    $local:path = [Environment]::GetEnvironmentVariable("Process", $scope).split(";")
    foreach ($scope in ("Machine","User")) {
        foreach ($var in [Environment]::GetEnvironmentVariable("Path", $scope).split(";")) {
            if (!($local:path -contains $var)) {
                $local:path += $var
            }
        }
    }
    $local:path = [system.String]::Join(";", $local:path)
    [Environment]::SetEnvironmentVariable("Path", $local:path, "Process")
}


# ============resetBinPath============
# Description: reset environement variable Path
function resetEnvPath() {
    $local:path = @()
    foreach ($scope in ("Machine","User")) {
        foreach ($var in [Environment]::GetEnvironmentVariable("Path", $scope).split(";")) {
            if (!($local:path -contains $var)) {
                $local:path += $var
            }
        }
    }
    $local:path = [system.String]::Join(";", $local:path)
    [Environment]::SetEnvironmentVariable("Path", $local:path, "Process")
}
    

# ============resetEnv============
# Description: reset environement variables
function resetEnv() {

    foreach ($key in [Environment]::GetEnvironmentVariables("Process").Keys) {
        [Environment]::SetEnvironmentVariable($key, $null, "Process")
    }
    foreach ($scope in ("Machine","User")) {
        $local:envvars = [Environment]::GetEnvironmentVariables($scope)
        foreach ($key in $local:envvars.Keys) {
            [Environment]::SetEnvironmentVariable($key, $local:envvars.Item($key), "Process")
        }
    }
    resetBinPath
}
