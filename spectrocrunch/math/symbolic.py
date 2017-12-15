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


class Function(object):

    @staticmethod
    def mergekwargs(kwargs1,kwargs2):
        kwargs = kwargs1.copy()
        kwargs.update(kwargs2)
        return kwargs
        
    def __init__(self,func,funcname=None,arglist=None,lowest=True,**kwargs):
        self._func = func
        self._funcname = funcname
        self._arglist = arglist
        self._lowest = lowest
        self._kwargs = kwargs

    def __call__(self,*args,**kwargs):
        _kwargs = self.mergekwargs(self._kwargs,kwargs)
        return self._func(*args,**_kwargs)

    def __str__(self):
        if self._funcname is None:
            name = self._func.__name__
        else:
            name = self._funcname
        if self._arglist is None:
            if self._lowest:
                arglist = "(...)"
            else:
                arglist = ""
        else:
            arglist = self._arglist
            
        return "{}{}".format(name,arglist)

    def __mul__(self,right):
        if isinstance(right,self.__class__):
            func = lambda *args,**kwargs: self._func(*args,**kwargs)*right._func(*args,**kwargs)
            name = "{}*{}".format(str(self),str(right))
        else:
            func = lambda *args,**kwargs: self._func(*args,**kwargs)*right
            name = "{}*{}".format(right,str(self))
        return self.__class__(func,name,lowest=False)
        
    def __rmul__(self,left):
        return self.__mul__(left)


