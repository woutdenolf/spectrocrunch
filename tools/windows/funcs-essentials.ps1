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
        netsh winhttp set proxy "http://proxy.esrf.fr:3128"
        (New-Object System.Net.WebClient).Proxy.Credentials = [System.Net.CredentialCache]::DefaultNetworkCredentials
    }
}


# ============download_file============
# Description: 
$webclient = New-Object System.Net.WebClient
function download_file([string]$url, [string]$output) {
	# Downloads a file if it doesn't already exist
	if (!(Test-Path $output -pathType leaf)){
        cprint "Downloading $url to $output ..."
		$webclient.DownloadFile($url, $output)
	}
}


# ============download_git_release============
# Description: 
function download_git_release([string]$user,[string]$project,[string]$extension,[AllowNull()][string]$output)
{
    $local:link = "https://api.github.com/repos/$user/$project/releases/latest"
    $response = Invoke-RestMethod -Method 'Get' -ContentType 'application/json' -Uri $local:link
    $filename = $null
    foreach ($asset in $response.assets) {
        if ($asset.name.endswith($extension)) {
            $filename = $asset.name
            break
        }
    }

    if ($filename -ne $null) {
        if ($output -ne $null -and $output -ne "") {
            $filename = $output
        }
        download_file $asset.browser_download_url $filename
    }

    return $filename
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