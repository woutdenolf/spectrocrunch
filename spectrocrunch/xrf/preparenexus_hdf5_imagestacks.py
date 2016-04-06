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

import numpy as np

def makeinfogroup(name,data):
    return {"name":name, "data":data}

def addinfogroup(fout,name,datadict):
    # Create info group when not there
    if "processing" not in fout:
        ginfo = fout.create_group("processing")
    else:
        ginfo = fout["processing"]

    # Add new info group
    name = "%d.%s"%(len(ginfo.keys())+1,name)
    newgroup = ginfo.create_group(name)
    for k in datadict:
        if isinstance(datadict[k],str):
            asciilist = [datadict[k].encode("ascii","ignore")]
            typ = "S%d"%len(asciilist[0])
            newgroup.create_dataset(k,(1,1),typ,asciilist)

        elif isinstance(datadict[k],list):
            arr = np.array(datadict[k])
            if arr.size==0:
                newgroup[k] = arr
            elif isinstance(arr[0],str) or isinstance(arr[0],np.unicode_):
                asciilist = [s.encode("ascii","ignore") for s in datadict[k]]
                typ = "S%d"%max([len(s) for s in asciilist])
                newgroup.create_dataset(k,(len(asciilist),1),typ,asciilist)
            else:
                newgroup[k] = arr
        else:
            newgroup[k] = datadict[k]

def addinfogroups(fout,infolist):
    for v in infolist:
        addinfogroup(fout,v["name"],v["data"])
    
def getinfogroups(fout):
    # Prepare list to receive info
    if "processing" not in fout:
        return []
    ginfo = fout["processing"]
    steps = ginfo.keys()
    if len(steps) == 0:
        return []
    ret = [None]*len(steps)

    # Get info
    for step in steps:
        tmp = step.split(".")
        ind = int(tmp[0])-1
        name = ".".join(tmp[1:])
        data = {k:ginfo[step][k][:] for k in ginfo[step]}
        ret[ind] = makeinfogroup(name,data)
    
    return ret

def copyaddinfogroup(fin,fout,name,datadict):
    #infolist = getinfogroups(fin)
    #infolist.append(makeinfogroup(name,datadict))
    #addinfogroups(fout,infolist)
    
    fin.copy(fin["processing"],fout)
    addinfogroup(fout,name,datadict)

def newgroups(fout,groups,extension):
    """Create new NXdata groups.

    Args:
        fout (hdf5.File): hdf5 file handle
        groups (list): group names without extension (str)
        extension (str): group name extension

    Returns:
        list: new NXdata groups (hdf5.Group)
    """

    if extension != "":
        add = "." + extension
    else:
        add = ""

    ngroups = len(groups)
    ret = [None]*ngroups
    for i in range(len(groups)):
        # Create the NXdata group
        nxdatagrp = fout.create_group(groups[i].split(".")[0] + add)
        nxdatagrp.attrs["NX_class"] = "NXdata"

        # Add dataset name
        nxdatagrp.attrs["signal"] = "data"
        ret[i] = nxdatagrp
    return ret

def newaxes(fout,axes,axesdata,extension):
    """Create new axes group.

    Args:
        fout (hdf5.File): hdf5 file handle
        axes (list): list of axes names without extension
        axesdata (list): new axes data
        extension (str): group name extension

    Returns:
        list: new axes names (dict)
    """

    if extension != "":
        add = "." + extension
    else:
        add = ""

    naxes = len(axes)
    ret = [None]*naxes

    for i in range(naxes):
        newname = axes[i]["fullname"].split(".")[0] + add
        if newname not in fout:
            fout[newname] = axesdata[i]
        ret[i] = {"fullname":newname,"name":axes[i]["name"]}

    return ret

def linkaxes(fout,axes,groups):
    """Link axes in NXdata groups

    Args:
        fout (hdf5.File): hdf5 file handle
        axes (list): list of axes names
        groups (list): list of NXdata groups
    """
    for nxdatagrp in groups:
        nxdatagrp.attrs["axes"] = [a["name"] for a in axes]
        for a in axes:
            if a["name"] in nxdatagrp:
                nxdatagrp[a["name"]].path = fout[a["fullname"]].name
            else:
                nxdatagrp[a["name"]] = fout[a["fullname"]]

