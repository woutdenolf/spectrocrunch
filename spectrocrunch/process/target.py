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


from abc import ABCMeta,abstractmethod
from future.utils import with_metaclass
from ..io import nxfs


class Target(with_metaclass(ABCMeta,object)):
    """task output
    """
    
    def __init__(self):
        self._counter = 0
        
    @property
    def name(self):
        num = self._counter0 + self._counter
        return self._fmt.format(num)
        
    @name.setter
    def name(self,value):
        # value == 'align':
        #   fmt = 'align.{:d}'
        #   counter0 = 1
        # value == 'align0002':
        #   fmt = 'align{:04d}'
        #   counter0 = 2
        self._fmt,self._counter0,nonumber = nxfs.targetname_fmt(value)
        if nonumber:
            self._counter0 = 1

    def increment(self):
        self._counter0 += 1
        
    @property
    @abstractmethod
    def exists(self):
        pass


class TargetFs(Target):
    """task output composed of a single filesystem node
    """

    def __init__(self,path):
        """
        Args:
            path(fs.Path)
        """
        super(TargetFs,self).__init__()
        self.path = path

    @property
    def path(self):
        return self.parent[self.name]
    
    @path.setter
    def path(self,value):
        self.parent = value.parent
        self.name = value.name
    
    @property
    def exists(self):
        return self.path.exists


class TargetFSArray(Target):
    """task output composed of multiple filesystem nodes
    """

    def __init__(self,parent,name,extensions):
        """
        Args:
            parent(fs.Path)
            name(str)
            extensions(list(str))
        """
        super(TargetFSArray,self).__init__()
        self.parent = parent
        self.name = name
        self.extensions = extensions

    def __iter__(self):
        name = self.name
        for ext in self.extensions:
            yield self.parent[name+ext]

    def __getitem__(self,ext):
        if ext in self.extensions:
            return self.parent[self.name+ext]
        else:
            raise IndexError('Invalid extension (should be {})'.format(self.extensions))
            
    @property
    def exists(self):
        return all(node.exists for node in self)
