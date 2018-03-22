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

from . import types
from ..common import instance 

import numpy as np

class StoichiometryBase(object):
    # TODO: refactor most of Mixture and Compound
        
    def setfraction(self,parts,values,fractype):
        if instance.isstring(parts):
            parts = [parts]
        values = instance.asarray(values)
    
        for p in parts:
            if p not in self.parts:
                raise RuntimeError("{} not in {}".format(p,self))
    
        # rebalance others
        w = self.fractions(fractype)
        w2 = dict(w)
        for p in parts:
            w2.pop(p)
        
        # update others
        if w2:
            v2 = np.asarray(w2.values())
            v2 *= (1-values.sum())/v2.sum()
            w.update((k,v) for k,v in zip(w2.keys(),v2))
        
        # update fractions
        w.update(zip(parts,values))
        self.change_fractions(w,fractype)
    
    def setmassfraction(self,comp,value):
        self.setfraction(comp,value,types.fraction.mass)

    def setmolefraction(self,comp,value):
        self.setfraction(comp,value,types.fraction.mole)
        
    def setvolumefraction(self,comp,value):
        self.setfraction(comp,value,types.fraction.volume)
        
