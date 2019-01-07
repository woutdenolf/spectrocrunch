#
# Essentials for installation and downloading
#


function require_build_essentials()
{
}


function require_web_essentials()
{
    require_web_access
}


function require_web_access()
{
    [Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12
    
    if ((Get-WmiObject Win32_ComputerSystem).Domain -eq "esrf.fr") {
        $global:http_proxy = "http://proxy.esrf.fr:3128"
        netsh winhttp set proxy $global:http_proxy
        $global:https_proxy = $global:http_proxy
        [Environment]::SetEnvironmentVariable("http_proxy", $global:http_proxy, "User")
        [Environment]::SetEnvironmentVariable("https_proxy", $global:http_proxy, "User")
        (New-Object System.Net.WebClient).Proxy.Credentials = [System.Net.CredentialCache]::DefaultNetworkCredentials
    }
}


# ============download_file============
# Description: 
$webclient = New-Object System.Net.WebClient
function download_file([string]$url, [AllowNull()][string]$output) {
    if ($output -eq $null -or $output -eq "") {
        $output = joinPath (Get-Location).Path ($url.split("/"))[-1]
    }
	# Downloads a file if it doesn't already exist
	if (!(Test-Path $output -pathType leaf)){
        cprint "Downloading $url to $output ..."
		$webclient.DownloadFile($url, $output)
    }
    return $output
}


# ============download_git_release============
# Description: download an asset that matches the extension pattern
function download_git_release([string]$user,[string]$project,[string]$extpattern,[AllowNull()][string]$output=$null,[string]$release="latest")
{
    $local:link = "https://api.github.com/repos/$user/$project/releases/$release"
    $response = Invoke-RestMethod -Method 'Get' -ContentType 'application/json' -Uri $local:link
    $filename = $null
    foreach ($asset in $response.assets) {
        $local:m = [regex]::match($asset.name,$local:extpattern)
        if ($local:m.Success) {
            $filename = joinPath (Get-Location).Path $asset.name
            break
        }
    }

    if ($filename -ne $null) {
        if ($output -ne $null -and $output -ne "") {
            $filename = $output
        }
        $filename = download_file $asset.browser_download_url $filename
    }

    return $filename
}


# ============yes/no============
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


# ============Unzip============
# Description: 
Add-Type -AssemblyName System.IO.Compression.FileSystem
function Unzip([string]$zipfile, [AllowNull()][string]$outpath)
{
    if ($outpath -eq $null -and $outpath -eq "") {
        $outpath = Get-Location
    }
    [System.IO.Compression.ZipFile]::ExtractToDirectory($zipfile, $outpath)
}


# ============Save-Env============
# Description: 
function Save-Env() {
    $tempFile = [IO.Path]::GetTempFileName()  
    Get-ChildItem env: | Export-CliXml $tempFile
    return $tempFile
}


# ============Restore-Env============
# Description: 
function Restore-Env([string]$tempFile) {
    Get-ChildItem env: | Foreach-Object { 
        set-item "env:\$($_.Name)" $null -force
    }
    Import-CliXml $tempFile | ForEach-Object {
        set-item "env:$($_.Name)" $_.Value -force
    }
    Remove-Item $tempFile -Force -ErrorAction SilentlyContinue
}


# ============Reset-Env============
# Description: reset environment variables
function Reset-Env([string]$scriptName) {
    $local:sep = "-=-cut-=-"
    $local:cmdLine = "echo $local:sep & set"
    if ([Environment]::Is64BitProcess) {
        $local:out = & $Env:SystemRoot\SysWOW64\cmd.exe /c $local:cmdLine | Out-String
    } else {
        $local:out = & $Env:SystemRoot\system32\cmd.exe /c $local:cmdLine | Out-String
    }
    if ($LASTEXITCODE -ne 0) {
        return
    }

    # Clear environment
    Get-ChildItem env: | Foreach-Object { 
        set-item "env:\$($_.Name)" $null -force
    }

    # Reset environment
    $local:arr = $local:out -split $local:sep
    if ($local:arr[0] -ne $null) {
        Write-Host $local:arr[0]
    }
    if ($local:arr[1] -ne $null -and $local:retcode -eq 0) {
        $local:arr[1].split([Environment]::NewLine) | select-string '^([^=]*)=(.*)$' | foreach-object {
            $varName = $_.Matches[0].Groups[1].Value
            $varValue = $_.Matches[0].Groups[2].Value
            set-item Env:$varName $varValue -force
            }
    }
}


# ============Show-Env============
# Description: 
function Show-Env() 
{
    Get-ChildItem env:
}


# ============Invoke-CmdScript============
# Description: run the batch script and inherit the environment
function Invoke-CmdScript([string]$scriptName) {
    $local:cmdLine = """$scriptName"" $args"
    $local:sep = "-=-cut-=-"
    $local:cmdLine = "$local:cmdLine & echo $local:sep & set"
    if ([Environment]::Is64BitProcess) {
        $local:out = & $Env:SystemRoot\SysWOW64\cmd.exe /c $local:cmdLine | Out-String
    } else {
        $local:out = & $Env:SystemRoot\system32\cmd.exe /c $local:cmdLine | Out-String
    }
    $local:retcode = $LASTEXITCODE
    $local:arr = $local:out -split $local:sep
    if ($local:arr[0] -ne $null) {
        Write-Host $local:arr[0]
    }
    if ($local:arr[1] -ne $null -and $local:retcode -eq 0) {
        $local:arr[1].split([Environment]::NewLine) | select-string '^([^=]*)=(.*)$' | foreach-object {
            $varName = $_.Matches[0].Groups[1].Value
            $varValue = $_.Matches[0].Groups[2].Value
            set-item Env:$varName $varValue -force
            }
    }
    return $local:retcode
}


# ============make-link============
# Description: 
function make-link($target,$link) {
    New-Item -ItemType SymbolicLink -Name $link -Target $target
}

# ============registry-value============
# Description: 
function registry-value($key,$property) {
    $local:keep = $ErrorActionPreference
    $ErrorActionPreference = "stop"
    try {
        return (Get-ItemProperty -Path $key -Name $property).$property
    } catch {
        return $null
    } finally {
        $ErrorActionPreference = $keep
    }
}