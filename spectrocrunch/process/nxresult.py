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

from ..utils.hashable import Hashable
from ..utils import instance
from ..io import nxfs
from . import axis

class Group(Hashable):
    """NXprocess subgroup name wrapper
    """
    
    def __init__(self,groupname):
        name = "counters"
        number = -1
        category = -1
            
        if instance.isnumber(groupname):
            name = "detector{:02d}".format(groupname)
            number = groupname
            category = 0
        elif instance.isstring(groupname):
            if groupname.startswith("xia"):
                groupname = groupname[3:]
            if groupname.startswith("detector"):
                groupname = groupname[8:]
                
            if groupname.isdigit():
                number = int(groupname)
                name = "detector{:02d}".format(number)
                category = 0
            elif groupname.startswith("S"):
                number = int(groupname[1:])
                name = "detectorS{:01d}".format(number)
                category = 1

        elif isinstance(groupname,self.__class__):
            name,number,category = groupname.name,groupname.number,groupname.category
        elif groupname:
            raise ValueError("Unexpected detector name {}".format(groupname))
        
        self.name = name
        self.number = number
        self.category = category

    @property
    def isdetector(self):
        # e.g. xia00, xiaS0
        return self.category>=0
        
    @property
    def issum(self):
        # e.g. xiaS0
        return self.category==1

    @property
    def issingledetector(self):
        # e.g. xia00
        return self.category==0

    @property
    def xialabel(self):
        if self.issingledetector:
            return 'xia{:02d}'.format(self.number)
        elif self.issum:
            return 'xiaS{:d}'.format(self.number)
        else:
            return None

    def __str__(self):
        return self.name

    def __repr__(self):
        return self.name

    def __eq__(self,other):
        if isinstance(other,self.__class__):
            return self.name==other.name
        else:
            return self.name==other
    
    def __ne__(self,other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(repr(self))


def regulargriddata(nxgroup):
    """
    Args:
        nxgroup(nxfs.Path): NXdata or NXprocess
    
    Returns:
        groups(dict): Group:list(nxfs.Path)
        axes(list(Axis)):
        stackdim(int):
    """
    axes = []
    groups = {}
        
    if nxgroup.nxclass=='NXdata':
        it = [nxgroup]
        stackdim = None
    elif nxgroup.nxclass=='NXprocess':
        progname = nxgroup['program'].read()
        if progname == nxfs.PROGRAM_NAME:
            stackdim = nxgroup.config.read().get('stackdim', None)
        else:
            raise ValueError('NXprocess from program "{}" is not known'.format(progname))
        it = nxgroup.results.iter_is_nxclass('NXdata')
    else:
        raise ValueError('{} should be an NXdata or NXprocess group'.format(nxgroup))

    for nxdata in it:
        if nxdata.islink:
            continue
        if stackdim is None:
            stackdim = nxdata.stackdim()
        group = Group(nxdata.name)
        if group in groups:
            raise RuntimeError('Group {} appears more than once'.format(group))
        axs = [axis.factory(values,name=name,title=attrs['title'],type='quantitative')
               for name,values,attrs in nxdata.axes]
        if axes:
            if axes!=axs:
                raise RuntimeError('{} has different axes'.format(nxdata))
        else:
            axes = axs
        groups[group] = list(nxdata.signals)

    return groups, axes, stackdim
