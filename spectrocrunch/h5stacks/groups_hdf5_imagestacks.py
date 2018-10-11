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

class Group(Hashable):

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
            elif groupname==name:
                pass
            elif groupname:
                raise ValueError("Unexpected detector name {}".format(groupname))
        elif isinstance(groupname,self.__class__):
            name,number,category = groupname.name,groupname.number,groupname.category
        elif groupname:
            raise ValueError("Unexpected detector name {}".format(groupname))
        
        self.name = name
        self.number = number
        self.category = category
            
    def _stringrepr(self):
        return self.name

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

    def __eq__(self,other):
        if isinstance(other,self.__class__):
            return self.name==other.name
        else:
            return self.name==other
    
    def __ne__(self,other):
        return not self.__eq__(other)
        
