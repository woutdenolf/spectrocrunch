# 
# Helper functions.
# 

function strreplace([string]$str,[string]$a,[string]$b)
{
    return $str.Replace("$a","$b")
}

function strtofloat([string]$str)
{
    return [float](strreplace $str "," ".")
}

function strtoint([string]$str)
{
    return [int](strreplace $str "," ".")
}
