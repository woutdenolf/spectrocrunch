# 
# This script will run Appveyor tests.
# 
# https://mnaoumov.wordpress.com/2015/01/11/execution-of-external-commands-in-powershell-done-right
#
# stdout is a string but stderr is an ErrorRecord,
# which raises NativeCommandError when written to stdout.
#
# Example: command that writes "hello" to stdout,
#          "world" to stderr and return with error code 1.
#
#   $ErrorActionPreference = "Continue"
#   $result = cmd /c "echo hello & echo world 1>&2 & exit 1" 2>&1
#
#  This returns [str, ErrorRecord] but would still raise an
#  exception if not for $ErrorActionPreference = "Continue".


function psexec([ScriptBlock] $ScriptBlock)
{
    $backupErrorActionPreference = $script:ErrorActionPreference
    $script:ErrorActionPreference = "Continue"
    try
    {
        # | % { ... }: foreach shorthand with $_ the loop variable
        & $ScriptBlock 2>&1 | % { "$_" }
    }
    finally
    {
        $script:ErrorActionPreference = $backupErrorActionPreference
    }

    if ($LastExitCode -ne 0) {$host.SetShouldExit($LastExitCode)}
}


function appveyor_unittest()
{
    cd $env:APPVEYOR_BUILD_FOLDER
    $PROJECTNAME = python setup.py name | Select-Object -Last 1

    cd $env:HOME
    psexec {python -m unittest -v "$PROJECTNAME.tests.test_all.main_test_suite"}
}

function appveyor_styletest()
{
    cd $env:APPVEYOR_BUILD_FOLDER
    psexec {python -m flake8 spectrocrunch}
}

function main()
{
    if ($env:APPVEYORRUN -eq "unit") {
        appveyor_unittest
    } elseif ($env:APPVEYORRUN -eq "style") {
        appveyor_styletest
    } else {
        Write-Host "No tests to run"
    }
}

main
