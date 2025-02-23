def search(names, name):
    lname = name.lower()
    # Contains (case insensitive)
    ret = [k for k in names if lname in k.lower()]
    if len(ret) > 1:
        # Equal (case insensitive)
        ret2 = [k for k in names if lname == k.lower()]
        if len(ret2) > 1:
            # Equal (case sensitive)
            ret3 = [k for k in names if name == k]
            if ret3:
                ret2 = ret3
        if ret2:
            ret = ret2
    return ret
