def search(names, name):
    name = name.lower()
    ret = [k for k in names if name in k.lower()]
    if len(ret) > 1:
        # Try exact match
        ret2 = [k for k in names if name == k.lower()]
        if ret2:
            ret = ret2
    return ret
