# -*- coding: utf-8 -*-
#
#   Copyright (C) 2017 European Synchrotron Radiation Facility, Grenoble, France
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


from ..common import listtools

from ..common.instance import isarray

import numpy as np

class discrete(object):
    
    def __init__(self,lines,intensities=None):
        """
        Args:
            lines(array(ureg.Quantity)): in keV, nm, ...
            intensities(array): line intensities
        """
        self.bnumber = not isarray(lines)

        if self.bnumber:
            lines = [lines]
    
        # Parse input
        if intensities is None:
            n = len(lines)
            intensities = [1./n]*n
        else:
            if isarray(intensities):
                self.bnumber = 0
            else:
                intensities = [intensities]

        # Sum equal lines
        lines,intensities = listtools.weightedsum(lines,intensities)
    
        # Sort by line
        self._lines,self._intensities = listtools.sort2lists(lines,intensities)
        
    def __add__(self, val):
        if isinstance(val,self.__class__):
            return self.__class__(self.lines + val.lines,\
                                  self.intensities + val.intensities)
        else:
            raise NotImplementedError
    
    @property
    def total(self):
        return sum(self._intensities)
    
    @property
    def size(self):
        return len(self._lines)
     
    @property   
    def ratios(self):
        if self.bnumber:
            return self._intensities[0]
        else:
            return np.asarray(self._intensities)/float(self.total)
    
    @property
    def lines(self):
        if self.bnumber:
            return self._lines[0]
        else:
            return self._lines
            
    @property
    def energies(self):
        if self.bnumber:
            return self._lines[0].to('keV','spectroscopy').magnitude
        else:
            return np.vertorize(lambda x:x.to('keV','spectroscopy').magnitude,otypes=[float])(self._lines)

        
