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

import h5py
from datetime import datetime
import numpy as np
import os
import pkg_resources

def hdf5pathparse(path):
    """
    Args:
        path(str):

    Returns:
        tuple: (path,name)
    """
    l = [x for x in path.split('/') if x]
    if len(l)==0:
        return '/',""
    elif len(l)==1:
        return '/',l[0]
    else:
        return '/'+'/'.join(l[:-1]),l[-1]

def timestamp():
    return datetime.now().isoformat()
    #return "T".join(str(datetime.now()).split())

class File(h5py.File):

    def __init__(self,filename,**kwargs):
        """
        r: readonly, file must exist
        r+: read/write, file must exist
        w: create, create when doesn't exist
        w-: create, fail if exists
        a: read/write when exists, create otherwise

        Raises:
            IOError: when file doesn't exist
        """
        h5py.File.__init__(self,filename,**kwargs)

        if "NX_class" not in self.attrs and self.mode!="r":
            self.attrs["NX_class"] = "NXroot"
            self.attrs["file_name"] = filename
            self.attrs["file_time"] = timestamp()
            self.attrs["HDF5_Version"] = h5py.version.hdf5_version
            self.attrs["h5py_version"] = h5py.version.version
            self.attrs["creator"] = "spectrocrunch"
            self.updated = True
        else:
            self.updated = False
            
    def markupdated(self):
        self.updated = True

    def close(self):
        if self.updated:
            self.attrs["file_update_time"] = timestamp()
        h5py.File.close(self)

def createlink(f,dest,linkdir,linkname,soft=True):
    """
    Args:
        f (h5py.File|str): hdf5 file
        dest (str): destination for the link
        linkdir (str): link directory
        linkname (str): link name
    """
    bclose = False
    if isinstance(f,h5py.File) or isinstance(f,h5py.Group):
        hdf5FileObject = f
    elif isinstance(f,str):
        hdf5FileObject = h5py.File(f)
        bclose = True
    else:
        raise ValueError("The hdf5 file must be either a string or an hdf5 file object.")

    if dest in hdf5FileObject:
        if soft:
            # Remove the link if it exists
            if linkdir in hdf5FileObject:
                if linkname in hdf5FileObject[linkdir]:
                    hdf5FileObject[linkdir].id.unlink(linkname)
            # Create the link
            hdf5FileObject[linkdir][linkname] = h5py.SoftLink(dest)
        else:
            b = True
            if linkdir in hdf5FileObject:
                if linkname in hdf5FileObject[linkdir]:
                    hdf5FileObject[linkdir][linkname].path = f[dest]
                    b = False
            if b:
                hdf5FileObject[linkdir][linkname] = f[dest]

    if bclose:
        hdf5FileObject.close()

def removesoftlinks(grp):
    """
    Args:
        grp(h5p.Group): group from which all links should be removed
    """
    for linkname in grp:
        grp.id.unlink(linkname)

def makeinfogroup(name,data):
    return {"name":name, "data":data}

def addinfogroup(fout,name,datadict):
    # Create info group when not there
    if "processing" not in fout:
        ginfo = newNXentry(fout,"processing")
    else:
        ginfo = fout["processing"]

    # Add new info group
    index = len(ginfo.keys())+1

    name = "%d.%s"%(index,name)
    newgroup = ginfo.create_group(name)
    newgroup.attrs["NX_class"] = "NXprocess"
    newgroup.attrs["program"] = "spectrocrunch"
    newgroup.attrs["version"] = pkg_resources.require("SpectroCrunch")[0].version
    newgroup.attrs["sequence_index"] = index
    newgroup.attrs["date"] = timestamp()

    for k in datadict:
        if (datadict[k] is None):
            continue    
        elif isinstance(datadict[k],str):
            asciilist = [datadict[k].encode("ascii","ignore")]
            typ = "S%d"%len(asciilist[0])
            newgroup.create_dataset(k,(1,1),typ,asciilist)
        elif isinstance(datadict[k],list) or isinstance(datadict[k],tuple):
            arr = np.array(datadict[k])
            arr[np.equal(arr,None)] = 0
            if arr.size==0:
                newgroup[k] = arr
            elif isinstance(arr.flat[0],str) or isinstance(arr.flat[0],np.unicode_):
                asciilist = [s.encode("ascii","ignore") for s in datadict[k]]
                typ = "S%d"%max([len(s) for s in asciilist])
                newgroup.create_dataset(k,(len(asciilist),1),typ,asciilist)
            elif isinstance(arr.flat[0],int):
                arr = np.array(arr,dtype=int)
            else:
                newgroup[k] = arr
        else:
            newgroup[k] = datadict[k]

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
        ret[ind] = {k:ginfo[step][k].value for k in ginfo[step]}
    
    return ret

def copyaddinfogroup(fin,fout,name,datadict):
    fin.copy(fin["processing"],fout)
    addinfogroup(fout,name,datadict)

def NXclassup(child,default=True,defaultup=True):
    parent = child.parent
    while "NX_class" not in parent.attrs:
        grandparent = parent.parent

        if isinstance(child,h5py.Dataset):
            parent.attrs["NX_class"] = "NXdata"
        else:
            # Parent is "NXentry" when grandparent is "NXroot",
            # otherwise parent is "NXsubentry".
            bentry = False
            if "NX_class" in grandparent.attrs:
                bentry = grandparent.attrs["NX_class"]=="NXroot"

            if bentry:
                parent.attrs["NX_class"] = "NXentry"
            else:
                parent.attrs["NX_class"] = "NXsubentry"

        child = parent
        parent = grandparent

def NXdefault(child,up=True):
    if isinstance(child,h5py.Dataset):
        path,name = hdf5pathparse(child.name) 
        p = child.parent
        p.attrs["signal"] = name
    else:
        p = child

    if not up:
        return

    # Update all "default" attributes of the parents
    while p.name != '/':
        path,name = hdf5pathparse(p.name)    
        p.parent.attrs["default"] = name
        p = p.parent

def defaultstack(f,nxdatagroup):
    """Make this NXdata group the default of the file

    Args:
        f (h5py.File|str): hdf5 file
        nxdatagroup (str): path of the default NXdata group
    """
    bclose = False
    if isinstance(f,h5py.File):
        hdf5FileObject = f
    elif isinstance(f,str):
        hdf5FileObject = h5py.File(f)
        bclose = True
    else:
        raise ValueError("The hdf5 file must be either a string or an hdf5 file object.")

    # Nexus default
    NXdefault(hdf5FileObject[nxdatagroup])

    # Pymca default: PyMcaGui/pymca/QHDF5StackWizard.py
    # It only expects 1 group (old txmwizard as well)
    #default = "_defaultstack"
    #if default in hdf5FileObject:
    #    entry = hdf5FileObject[default]
    #else:
    #    entry = newNXentry(hdf5FileObject,default)
    #entry.attrs["description"] = "For viewers who don't implement the new default data convention."
    #removesoftlinks(entry)
    #path, name = hdf5pathparse(nxdatagroup)
    #createlink(hdf5FileObject,nxdatagroup,entry.name,name)

    # Done
    if bclose:
        hdf5FileObject.close()

def newNXdata(parent,name,extension):
    """Create new NXdata groups.

    Args:
        parent (h5py.File or h5py.Group): hdf5 file handle or group
        name (str): group name without extension (str)
        extension (str): group name extension

    Returns:
        hdf5.Group: new NXdata groups
    """

    if extension != "":
        add = "." + extension
    else:
        add = ""

    # Create the NXdata group
    nameext = name.split(".")[0] + add
    nxdatagrp = parent.create_group(nameext)
    nxdatagrp.attrs["NX_class"] = "NXdata"
    nxdatagrp.attrs["signal"] = "data"

    # Entry/subentry all the way up till the root
    NXclassup(nxdatagrp)

    # Default group
    NXdefault(nxdatagrp)

    return nxdatagrp

def newNXentry(parent,name):
    """Create new NXentry group.

    Args:
        parent (h5py.File or h5py.Group): hdf5 file handle or group
        group (str): group name without extension (str)

    Returns:
        hdf5.Group: new NXentry group
    """

    grp = parent.create_group(name)
    grp.attrs["NX_class"] = "NXentry"
    if "NX_class" in parent.attrs:
        if parent.attrs["NX_class"] == "NXentry":
            grp.attrs["NX_class"] = "NXsubentry"
    return grp

def createNXdataSignal(nxdatagrp,**kwargs):
    """Create the signal dataset in a NXdata group

    Args:
        nxdatagrp(h5py.Group): NXdata group
    Returns:
        h5py.Dataset
    """

    if "data" in kwargs:
        if isinstance(kwargs["data"],h5py.Dataset):
            if kwargs["data"].file == nxdatagrp.file:
                # Dataset and NXdata group in the same file: create link
                createlink(nxdatagrp,kwargs["data"].name,nxdatagrp.name,nxdatagrp.attrs["signal"])
                return
            else:
                pass
                #kwargs["data"] = kwargs["data"].value

    # Create dataset
    dset = nxdatagrp.create_dataset(nxdatagrp.attrs["signal"],**kwargs)
    dset.attrs["signal"] = 1
    return dset

def createaxes(fout,axes):
    """Create new axes group

    Args:
        fout (h5py.File): hdf5 file handle
        axes (list(dict)): [{"name":"name1","data":np.array},
        {"name":"name2","data":np.array},
        {"name":"name3","data":np.array}]

    Returns:
        list(dict): [{"name":"name1","fullname":"/axes/name1/data"},
        {"name":"name2","fullname":"/axes/name2/data"},
        {"name":"name3","fullname":"/axes/name3/data"}]
    """
    if "axes" in fout:
        grp = fout["axes"]
    else:
        grp = newNXentry(fout,"axes")
    naxes = len(axes)
    ret = [None]*naxes
    for i in range(naxes):
        grp2 = newNXdata(grp,axes[i]["name"],"")
        dset = grp2.create_dataset(grp2.attrs["signal"],data=axes[i]["data"])
        #dset.attrs["axis"] = i+1 # THIS IS DEPRICATED
        ret[i] = {"name":axes[i]["name"],"fullname":dset.name}
    return ret

def newaxes(fout,axes,axesdata,extension):
    """Create new axes group if it doesn't exist yet.

    Args:
        fout (h5py.File): hdf5 file handle
        axes (list): list of axes names without extension
        axesdata (list): new axes data
        extension (str): group name extension

    Returns:
        list(dict): [{"name":"name1","fullname":"/axes/name1/data"},
                     {"name":"name2","fullname":"/axes/name2/data"},
                     {"name":"name3","fullname":"/axes/name3/data"}]
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
            if isinstance(axesdata[i],h5py.Dataset):
                path,name = hdf5pathparse(newname)
                createlink(fout,axesdata[i].name,path,name,soft=True)
            else:
                fout[newname] = axesdata[i]
                #fout[newname].attrs["axis"] = i+1 # THIS IS DEPRICATED
                NXclassup(fout[newname])
                NXdefault(fout[newname],up=False)
        ret[i] = {"fullname":newname,"name":axes[i]["name"]}

    return ret

def linkaxes(fout,axes,groups):
    """Link axes in NXdata groups

    Args:
        fout (h5py.File): hdf5 file handle
        axes (list(str)): list of axes names
        groups (list(str)): list of NXdata groups
    """
    for nxdatagrp in groups:
        axesstr = ':'.join([a["name"] for a in axes])
        nxdatagrp.attrs["axes"] = axesstr
        nxdatagrp[nxdatagrp.attrs["signal"]].attrs["axes"] = axesstr
        for i in range(len(axes)):
            a = axes[i]
            createlink(fout,a["fullname"],nxdatagrp.name,a["name"],soft=True)
            nxdatagrp.attrs[a["name"]+"_indices"]=i

def parse_NXdata(grp):
    axesnames = grp.attrs["axes"].split(':')
    data = grp[grp.attrs["signal"]]
    axes = [grp[a] for a in axesnames]
    return data,axes,axesnames
    
