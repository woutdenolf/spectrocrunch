#
# Git
#

. $PSScriptRoot\funcs.ps1
. $PSScriptRoot\funcs-essentials.ps1
. $PSScriptRoot\funcs-install.ps1

function require_git()
{
    cprint "Checking git ..."

    if (cmdexists git) {
        cprint "Git is installed"
        return
    }

    if ((install_arch) -eq 64) {
        $local:extension = "64-bit.exe"
        $local:affix = ""
        $local:name = "Git (64-bit)"
    } else {
        $local:extension = "32-bit.exe"
        $local:affix = "-32"
        $local:name = "Git (32-bit)"
    }

    cprint "Download git ..."
    if (!(dryrun)) {
        require_web_access
        $local:filename = download_git_release "git-for-windows" "git" $local:extension
    }

    cprint $local:filename
    $local:path = joinPath (project_prefix) "Git$affix"
    cprint "Installing $name in $path ..."
    if (!(dryrun) -and $local:filename -ne $null) {
        $arguments = @()
        $arguments += "/silent"
        $arguments += "/dir=""$path"""
        install_exe $filename $arguments
        prependBinPath (joinPath $path "bin")

        if (cmdexists git) {
            cprint "Git is installed."
        } else {
            cerror "Git is not installed."
        }
    }
}