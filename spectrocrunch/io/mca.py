# -*- coding: utf-8 -*-
#
#   Copyright (C) 2018 European Synchrotron Radiation Facility, Grenoble, France
#
#   Principal author:   Wout De Nolf (wout.de_nolf@esrf.eu)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import time
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
    
def save(mca,filename,mode="w",zero=0,gain=1):
    with open(filename, mode) as ffile:
        ffile.write("#F %s\n" % filename)
        ffile.write("#D %s\n" % (time.ctime(time.time())))
        ffile.write("\n")
        ffile.write("#S 1 %s\n" % "spectrocrunch")
        ffile.write("#D %s\n" % (time.ctime(time.time())))
        ffile.write("#@MCA %16C\n")
        ffile.write("#@CHANN %d %d %d 1\n" % (len(mca), 0, len(mca)-1))
        ffile.write("#@CALIB %.7g %.7g %.7g\n" % (zero,gain,0))
        ffile.write(array2SpecMca(mca))
        ffile.write("\n")

def read(filename):
    f = spec(filename)
    return f.getdata2(1,[])

        
