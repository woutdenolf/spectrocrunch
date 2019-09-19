# -*- coding: utf-8 -*-

import time
import os
from .utils import mkdir
from .spec import spec


def array2SpecMca(data):
    """ Write a python array into a Spec array.
        Return the string containing the Spec array
    """
    tmpstr = "@A "
    length = len(data)
    for idx in range(0, length, 16):
        if idx+15 < length:
            for i in range(0, 16):
                tmpstr += "%.8g " % data[idx+i]
            if idx+16 != length:
                tmpstr += "\\"
        else:
            for i in range(idx, length):
                tmpstr += "%.8g " % data[i]
        tmpstr += "\n"
    return tmpstr


def save(mca, filename, mode="w", zero=0, gain=1):
    mkdir(os.path.dirname(filename))
    with open(filename, mode) as ffile:
        ffile.write("#F %s\n" % filename)
        ffile.write("#D %s\n" % (time.ctime(time.time())))
        ffile.write("\n")
        ffile.write("#S 1 %s\n" % "spectrocrunch")
        ffile.write("#D %s\n" % (time.ctime(time.time())))
        ffile.write("#@MCA %16C\n")
        ffile.write("#@CHANN %d %d %d 1\n" % (len(mca), 0, len(mca)-1))
        ffile.write("#@CALIB %.7g %.7g %.7g\n" % (zero, gain, 0))
        ffile.write(array2SpecMca(mca))
        ffile.write("\n")


def read(filename):
    f = spec(filename)
    return f.getdata2(1, [])
