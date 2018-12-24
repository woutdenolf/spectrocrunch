#
# Microsoft Visual Studio compiler
#

function msc_versions()
{
    $versions = @{}

    $tmp = @{}
    $tmp["version"] = "6.0"
    $tmp["name"] = "Visual Studio 6.0"
    $versions[1200] = $tmp

    $tmp = @{}
    $tmp["version"] = "7.0"
    $tmp["name"] = "Visual Studio .NET 2002"
    $versions[1300] = $tmp

    $tmp = @{}
    $tmp["version"] = "7.1"
    $tmp["name"] = "Visual Studio .NET 2003"
    $versions[1310] = $tmp

    $tmp = @{}
    $tmp["version"] = "8.0"
    $tmp["name"] = "Visual Studio 2005"
    $versions[1400] = $tmp

    $tmp = @{}
    $tmp["version"] = "9.0"
    $tmp["name"] = "Visual Studio 2008"
    $versions[1500] = $tmp

    $tmp = @{}
    $tmp["version"] = "10.0"
    $tmp["name"] = "Visual Studio 2010"
    $versions[1600] = $tmp

    $tmp = @{}
    $tmp["version"] = "11.0"
    $tmp["name"] = "Visual Studio 2012"
    $versions[1700] = $tmp

    $tmp = @{}
    $tmp["version"] = "12.0"
    $tmp["name"] = "Visual Studio 2013"
    $versions[1800] = $tmp

    $tmp = @{}
    $tmp["version"] = "14.0"
    $tmp["name"] = "Visual Studio 2015"
    $tmp["vcvarsall"] = "C:\Program Files (x86)\Microsoft Visual Studio\2017\BuildTools\VC\Auxiliary\Build\vcvarsall.bat"
    $tmp[32] = "x86 -vcvars_ver=14.0"
    $tmp[64] = "x64 -vcvars_ver=14.0"
    $versions[1900] = $tmp

    $tmp = @{}
    $tmp["version"] = "15.0"
    $tmp["name"] = "Visual Studio 2017"
    $tmp["vcvarsall"] = "C:\Program Files (x86)\Microsoft Visual Studio\2017\BuildTools\VC\Auxiliary\Build\vcvarsall.bat"
    $tmp[32] = "x86 -vcvars_ver=15.0"
    $tmp[64] = "x64 -vcvars_ver=15.0"
    $versions[1910] = $tmp

    $tmp = @{}
    $tmp["version"] = "15.3"
    $tmp["name"] = "Visual Studio 2017"
    $tmp["vcvarsall"] = "C:\Program Files (x86)\Microsoft Visual Studio\2017\BuildTools\VC\Auxiliary\Build\vcvarsall.bat"
    $tmp[32] = "x86 -vcvars_ver=15.3"
    $tmp[64] = "x64 -vcvars_ver=15.3"
    $versions[1911] = $tmp

    $tmp = @{}
    $tmp["version"] = "15.5"
    $tmp["name"] = "Visual Studio 2017"
    $tmp["vcvarsall"] = "C:\Program Files (x86)\Microsoft Visual Studio\2017\BuildTools\VC\Auxiliary\Build\vcvarsall.bat"
    $tmp[32] = "x86 -vcvars_ver=15.5"
    $tmp[64] = "x64 -vcvars_ver=15.5"
    $versions[1912] = $tmp

    $tmp = @{}
    $tmp["version"] = "15.6"
    $tmp["name"] = "Visual Studio 2017"
    $tmp["vcvarsall"] = "C:\Program Files (x86)\Microsoft Visual Studio\2017\BuildTools\VC\Auxiliary\Build\vcvarsall.bat"
    $tmp[32] = "x86 -vcvars_ver=15.6"
    $tmp[64] = "x64 -vcvars_ver=15.6"
    $versions[1913] = $tmp

    $tmp = @{}
    $tmp["version"] = "15.7"
    $tmp["name"] = "Visual Studio 2017"
    $tmp["vcvarsall"] = "C:\Program Files (x86)\Microsoft Visual Studio\2017\BuildTools\VC\Auxiliary\Build\vcvarsall.bat"
    $tmp[32] = "x86 -vcvars_ver=15.7"
    $tmp[64] = "x64 -vcvars_ver=15.7"
    $versions[1914] = $tmp

    $tmp = @{}
    $tmp["version"] = "15.8"
    $tmp["name"] = "Visual Studio 2017"
    $tmp["vcvarsall"] = "C:\Program Files (x86)\Microsoft Visual Studio\2017\BuildTools\VC\Auxiliary\Build\vcvarsall.bat"
    $tmp[32] = "x86 -vcvars_ver=15.8"
    $tmp[64] = "x64 -vcvars_ver=15.8"
    $versions[1915] = $tmp

    $tmp = @{}
    $tmp["version"] = "15.9"
    $tmp["name"] = "Visual Studio 2017"
    $tmp["vcvarsall"] = "C:\Program Files (x86)\Microsoft Visual Studio\2017\BuildTools\VC\Auxiliary\Build\vcvarsall.bat"
    $tmp[32] = "x86 -vcvars_ver=15.9"
    $tmp[64] = "x64 -vcvars_ver=15.9"
    $versions[1916] = $tmp

    return $versions
}

function msc_info([int]$msc_ver)
{
    return $(msc_versions)[$msc_ver]
}
