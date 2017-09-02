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

from ..common.classfactory import FactoryMeta
from future.utils import with_metaclass
import numpy as np

class SimulClass(object):
    
    @staticmethod
    def raiseabstract():
        raise NotImplementedError("SimulClass is an abstract class and shouldn't be instantiated.")  
        
    @classmethod
    def required(cls,arg,strarg):
        if arg is None:
            raise RuntimeError("{} not defined for {}".format(strarg.capitalize(),cls.__name__))
    
    @staticmethod
    def defined(func,var,strvar):
        if var is None:
            raise RuntimeError("{} not defined for {}".format(strvar,func.__name__))

    @staticmethod
    def broadcastold(N,energy):
        """
        Args:
            N(num|array): incomming number of photons with uncertainties
            energy(num|array): associated energies
            
        Returns:
            unumpy.uarray: len(energy) x len(N)
        """
        nen = np.asarray(energy).size
        N = np.asarray(N)
        if N.ndim==2:
            return np.broadcast_to(N,[nen,N.shape[1]])
        return np.broadcast_to(N,[nen,N.size])

def with_simulmetaclass(bases=None):
    if bases is None:
        return with_metaclass(FactoryMeta,SimulClass)
    else:
        return with_metaclass(FactoryMeta,SimulClass,*bases)
