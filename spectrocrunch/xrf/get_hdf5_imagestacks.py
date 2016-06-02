# -*- coding: utf-8 -*-
#
#   Copyright (C) 2015 European Synchrotron Radiation Facility, Grenoble, France
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

import spectrocrunch.io.nexus as nexus

def get_hdf5_imagestacks(h5file,datagroupnames):
    f = nexus.File(h5file,mode='r')

    # Get axes
    axesdict = {}
    for k in f["axes"].keys():
        if "." in k:
            continue
        name = str(k)
        signal = f["axes"][k].attrs["signal"]
        axesdict[name] = {"fullname":f["axes"][k][signal].name,"name":name}
    axes = None

    # Get data groups
    groups = [k for k in f.keys() if k in datagroupnames]
    stacks = {}
    for grp in groups:
        stacks[grp] = {k:f[grp][k].name for k in f[grp].keys() if "." not in k}
        if axes is None:
            names = f[stacks[grp].values()[0]].attrs["axes"].split(':')
            axes = [axesdict[name] for name in names]

    f.close()

    return stacks, axes
