#!/bin/bash
# 
# Helper functions.
# 


function strlen()
{
    echo ${#1}
}


function substr()
{
    local a=${2}
    local b=${3}
    local n=$(strlen ${1})
    if [[ ${a} -lt 0 ]]; then
        a=$((a+n))
    fi
    if [[ ${b} -lt 0 ]]; then
        b=$((b+n))
    fi
    n=$((b-a+1))
    expr "${1:${a}:${n}}"
}


function nchar_occurence()
{
    echo "${1}" | awk -F"${2}" '{print NF-1}'
}

function strreplace()
{
    echo -e ${1//${2}/${3}}
}



