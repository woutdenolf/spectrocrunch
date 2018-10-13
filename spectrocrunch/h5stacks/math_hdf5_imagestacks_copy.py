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

from ..io import nexus

def evaluate_entire(operation,fin,varargs,retstacks):
    """Copy image stacks
    
    Args:
        operation(dict):
        fin(nexus.File):
        varargs(dict):
        fixedargs(dict):
        retstacks(list(h5py.Group)):
    """
    
    for i in range(len(varargs)):
        grp = fin[varargs[i]]
        data = grp[grp.attrs["signal"]]
        dset = nexus.createNXdataSignal(retstacks[i],data=data,chunks = True)

def evaluate_sliced(operation,fin,varargs,retstacks):
    """Copy image stacks, slice by slice (is this needed?)
    
    Args:
        operation(dict)
        fin(nexus.File)
        varargs(dict)
        fixedargs(dict)
        retstacks(list(h5py.Group))
    """
    
    stackdim = operation["stackdim"]

    for i in range(len(varargs)):
        grp = fin[varargs[i]]
        data = grp[grp.attrs["signal"]]
        dset = nexus.createNXdataSignal(retstacks[i],shape=data.shape,dtype=data.dtype,chunks = True)

        nslice = data.shape[stackdim]
        for j in range(nslice):
            if stackdim==0:
                dset[j,...] = data[j,...]
            elif stackdim==1:
                dset[:,j,:] = data[:,j,:]
            else:
                dset[...,j] = data[...,j]

def evaluate(operation,fin,varargs,retstacks):
    if "sliced" in operation:
        if operation["sliced"]:
            evaluate_sliced(operation,fin,varargs,retstacks)
            return
    evaluate_entire(operation,fin,varargs,retstacks)

