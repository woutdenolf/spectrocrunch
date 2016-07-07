# -*- coding: utf-8 -*-
#
#   Copyright (C) 2016 European Synchrotron Radiation Facility, Grenoble, France
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

import os
import spectrocrunch.io.nexus as nexus

def hdf5base(file_in):
    tmp,ext = os.path.splitext(file_in)
    base,_ = os.path.splitext(tmp)
    return base,ext

def isgroup(hdf5path,name):
    """Check whether the name corresponds with the hdf5 path
    """
    return hdf5path.split(".")[0].endswith(name)

def selectgroups(hdf5paths,names):
    """Select hdf5 paths which correspond to the names
    """
    return [f for f in hdf5paths if any(isgroup(f,s) for s in names)]

def selectnotgroups(hdf5paths,names):
    """Select hdf5 paths which do not correspond to the names
    """
    return [f for f in hdf5paths if not any(isgroup(f,s) for s in names)]

def getnames(hdf5paths):
    """Get name from hdf5 path (/detector0/Si-K.align -> Si-K)
    """
    return [s.rstrip("/").split("/")[-1].split(".")[0] for s in hdf5paths]

def defaultstack(f,stacks,default):
    """Make the default group in a nexus file
    """
    if default is None:
        reference = []
    else:
        reference = [s for s in stacks if isgroup(s,default)]
    if len(reference)!=1:
        return
    nexus.defaultstack(f,reference[0])

def flattenstacks(stacks):
    """ {'detector0': {'xmap_x2': '/detector0/xmap_x2',...}, ...} -> ['/detector0/xmap_x2',...]
    """
    ret = []
    for k in stacks:
        ret += stacks[k].values()
    return ret

