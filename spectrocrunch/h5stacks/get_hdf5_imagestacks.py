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

import ..io.nexus as nexus

import re

def get_hdf5_imagestacks(h5file,datagroupnames):
    """
    Args:
        h5file(str): hdf5 file
        datagroupnames(list(str)): NXentries containing NXdata image stacks

    Returns:
        tuple

        The first element contains the image stack:
            stacks = {"counters":{"name1":nxdatapath1,"name2":nxdatapath2,...},
                      "det0":{"name3":nxdatapath3,"name4":nxdatapath4,...},
                      "det1":{"name3":nxdatapath5,"name4":nxdatapath6,...},...}

        The second element is a list with three elements which contains
        the axis values of the stack:
            stackaxes = [{"name":"name1","fullname":"path1"},
                         {"name":"name2","fullname":"path2"},
                         {"name":"name3","fullname"":"path3"}]
    """

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
    groups = [k for k in f.keys() if any(re.match(pattern,k) is not None for pattern in datagroupnames)]
    stacks = {}
    for grp in groups:
        stacks[grp] = {k:f[grp][k].name for k in f[grp].keys() if "." not in k}
        if axes is None and len(stacks[grp])!=0:
            names = f[stacks[grp].values()[0]].attrs["axes"].split(':')
            axes = [axesdict[name] for name in names]

    f.close()

    return stacks, axes
