#
# powershell package manager
#

. $PSScriptRoot\funcs.ps1
. $PSScriptRoot\funcs-essentials.ps1

function require_ps_packagemgr()
{
    if ((cmdexists Get-VSSetupInstance)) {
        cprint "VSSetup is installed"
        return
    }

    require_web_access

    if ((cmdexists Install-Module)) {
        Install-Module VSSetup -Scope CurrentUser
    } else {
        $local:filename = download_git_release "Microsoft" "vssetup.powershell" ".zip"
        if ($local:filename -ne $null) {
            Unzip $local:filename "$([Environment]::GetFolderPath("MyDocuments"))\WindowsPowerShell\Modules\VSSetup"
        }
    }

    if ((cmdexists Get-VSSetupInstance)) {
        cprint "VSSetup is installed"
    } else {
        cerror "VSSetup is not installed"
    }
}
