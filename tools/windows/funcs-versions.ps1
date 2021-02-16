# ============require_new_version============
# Description: a new version is required when current < required
function require_new_version([AllowNull()][string]$currentversion,[AllowNull()][string]$requiredversion)
{
    if ($currentversion -eq $null -or $currentversion -eq "") {
        # not version not exists
        return $true
    }

    if ($requiredversion -eq $null -or $requiredversion -eq "") {
        # no specific version required
        return $false
    }

    $local:currentv = $currentversion.split(".")
    $local:requiredv = $requiredversion.split(".")
    $local:ncurrentv = $local:currentv.Length
    $local:nrequiredv = $local:requiredv.Length
    $local:n = [math]::min($local:ncurrentv,$local:nrequiredv)
    foreach ($i in $(0..($local:n-1))) {
        if ([int]$local:currentv[$i] -lt [int]$local:requiredv[$i])
        {
            return $true
        }
    }
    return $false
}

# ============extract_version============
# Description:
function extract_version([string]$s)
{
    $local:m = [regex]::match($s,"[\d\.]+[\d]")
    if ($local:m.Success) {
        return $local:m.Groups[0].Value
    }
}

# ============require_new_version_strict============
# Description: a new version is required when current != required (up to a common depth)
function require_new_version_strict([AllowNull()][string]$currentversion,[AllowNull()][string]$requiredversion)
{
    if ($currentversion -eq $null -or $currentversion -eq "") {
        # no version exists
        return $true
    }

    if ($requiredversion -eq $null -or $requiredversion -eq "") {
        # no specific version required
        return $false
    }

    $local:currentv = $currentversion.split(".")
    $local:requiredv = $requiredversion.split(".")
    $local:n = [math]::min($local:currentv.Length,$local:requiredv.Length)-1
    $local:currentv = [string]::Join("",$local:currentv[0..$local:n])
    $local:requiredv = [string]::Join("",$local:requiredv[0..$local:n])
    return $local:currentv -ne $local:requiredv
}


# ============get_local_version============
# Description: returns a local version when it equal or newer than the requested version
function get_local_version([AllowNull()][string]$requiredv)
{
    foreach ($path in Get-ChildItem) {
        if ($path.Attributes -ne "Directory") {
            $local:version = extract_version $path.Name
            if ($local:version -ne $null) {
                if (!(require_new_version $local:version $requiredv)) {
                    return ($local:version,$path.FullName)
                }
            }
        }
    }
}


# ============get_local_version_strict============
# Description: returns a local version when it matches the requested version (up to a common depth)
function get_local_version_strict([AllowNull()][string]$requiredv)
{
    foreach ($path in Get-ChildItem) {
        if ($path.Attributes -ne "Directory") {
            $local:version = extract_version $path.Name
            if ($local:version -ne $null) {
                if (!(require_new_version_strict $local:version $requiredv)) {
                    return ($local:version,$path.FullName)
                }
            }
        }
    }
}


# ============version_intarray============
# Description: 
function version_intarray([string]$version,[AllowNull()][int]$n)
{
    $local:tmp = $version.split('.')
    if ($n -ne $null) {
        if (($local:tmp.Length) -lt $n) {
            $local:tmp += @(0)*($n-($local:tmp.Length))
        }
        if (($local:tmp.Length) -gt $n) {
            $local:tmp = $local:tmp[0..($n-1)]
        }
    }
    return $local:tmp | % {iex $_}
}


# ============version_pattern============
# Description: 
function version_pattern([string]$version,[int]$n)
{
    if ($local:version -eq $null -or $local:version -eq "") {
        $local:pattern = @("[\d]+")*$n
    } else {
        $local:pattern = $local:version.split('.')
        if ($n -gt 0) {
            if ($local:pattern.Length -lt $n) {
                $local:pattern += @("[\d]+")*($n-$local:pattern.Length)
            } else {
                $local:pattern = $local:pattern[0..($n-1)]
            }
        }
    }
    if ($local:pattern[-1] -eq "[\d]+") {
        $local:pattern = [string]::Join("\.",$local:pattern)
        return "($local:pattern)"
    } else {
        $local:pattern = [string]::Join("\.",$local:pattern)
        return "($local:pattern)[^\d]"
    }
}


# ============version_pattern============
# Description: 
function version_intarray_higher($oldarr,$newarr)
{
    $local:n = $oldarr.Length
    foreach ($i in $(0..($local:n-1))) {
        if ($newarr[$i] -gt $oldarr[$i]) {
            return $true
        } elseif ($newarr[$i] -lt $oldarr[$i]) {
            return $false
        }
    }
    return $false
}


# ============git_release============
# Description: get release that matches version and extension pattern
function git_release([string]$user,[string]$project,[string]$extpattern,[AllowNull()][string]$version)
{
    if ($version -eq $null -or $version -eq "")
    {
        return $null
    }

    $local:n = 3
    $local:pattern = version_pattern $local:version $local:n
    $local:mversion = @(0)*$local:n
    $local:id = $null
    $local:link = "https://api.github.com/repos/$user/$project/releases"
    $response = Invoke-RestMethod -Method 'Get' -ContentType 'application/json' -Uri $local:link
    foreach ($release in $response) {
        $local:m = [regex]::match($release.tag_name,$local:pattern)
        if ($local:m.Success) {
            $local:hasext = $false
            foreach ($asset in $release.assets) {
                $local:hasext = $local:m -or (([regex]::match($asset.name,$extpattern)).Success)
            }
            if ($local:hasext) {
                $local:tmp = version_intarray $m.Groups[1].Value $local:n
                if ((version_intarray_higher $local:mversion $local:tmp)) {
                    $local:mversion = $local:tmp
                    $local:id = $release.id
                }
            }
        }
    }
    return $local:id
}


# ============get_version_link============
# Description: get release that matches version and extension pattern
# Avoid when possible (limited number of calls)
function get_version_link([string]$url,[string]$pattern1,[string]$pattern2,[int]$n=3) {
    $local:mversion = @(0)*$n
    $local:mlink = $null
    $local:version = $null
    $local:content1 = Invoke-WebRequest -Uri $url
    foreach ($link1 in $local:content1.Links) {
        $local:m = [regex]::match($link1.href,$local:pattern1)
        if ($local:m.Success) {
            $local:tmp = version_intarray $m.Groups[1].Value $n
            if ((version_intarray_higher $local:mversion $local:tmp)) {
                $local:content2 = Invoke-WebRequest -Uri "$url/$($link1.href)"
                foreach ($link2 in $local:content2.Links) {
                    $local:m = [regex]::match($link2.href,$local:pattern2)
                    if ($local:m.Success) {
                        $local:tmp = version_intarray $m.Groups[1].Value $n
                        if ((version_intarray_higher $local:mversion $local:tmp)) {
                            $local:mversion = $local:tmp
                            $local:mlink = "$url/$($link1.href)/$($link2.href)"
                            $local:version = $m.Groups[1].Value
                        }
                    }
                }
            }
        }
    }
    return ($local:version,$local:mlink)
}