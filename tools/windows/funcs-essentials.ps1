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
