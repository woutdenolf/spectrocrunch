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

from __future__ import division

import numbers

class linop(object):

    def __init__(self,m,b):
        self.m = m
        self.b = b

    def __call__(self,x):
        return self.m*x + self.b

    def __mul__(self,rother):
        # first self, then rother
        if isinstance(rother,self.__class__): 
            return self.__class__(self.m*rother.m,self.b*rother.m+rother.b)
        else:
            raise NotImplementedError

    def __div__(self,rother):
        # first self, then rother
        if isinstance(rother,self.__class__):
            return self.__class__(self.m/rother.m,self.b/rother.m-rother.b/rother.m)
        else:
            raise NotImplementedError

    def __rmul__(self,lother):
        # first lother, then self
        if isinstance(lother,self.__class__): 
            return lother.__mul__(self)
        else:
            raise NotImplementedError

    def __rdiv__(self,lother):
        # first lother, then self
        if isinstance(lother,self.__class__): 
            return lother.__div__(self)
        else:
            raise NotImplementedError

    def __pow__(self,p):
        if isinstance(p,numbers.Integral):
            ret = self.__class__(self.m,self.b)
            for i in range(1,abs(p)):
                ret *= self.__class__(self.m,self.b)
            if p<0:
                ret = self.__class__(1/ret.m,-ret.b/ret.m)
            return ret
        else:
            raise NotImplementedError

    def __eq__(self,other):
        return self.m==other.m and self.b==other.b

    def __ne__(self,other):
        return self.m!=other.m and self.b!=other.b

    def __str__(self):
        return "y = {} x{:+}".format(self.m,self.b)

    def __repr__(self):
        return "{}({},{})".format(self.__class__,self.m,self.b)
