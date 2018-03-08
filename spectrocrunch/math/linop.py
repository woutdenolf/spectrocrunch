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
import numpy as np
import copy
import operator

from ..common import units

# TODO: do this in sympy


class Operator(object):
    """An operator with an operator on the right
    """
    
    def __init__(self):
        self._opright = None # equivalent to the identity

    @property
    def _op_name(self):
        return self.__class__.__name__
    
    @property
    def _op_args(self):
        return ""
    
    def __str__(self):
        args = self._op_args
    
        if self._opright is None:
            var = "x"
        else:
            var = str(self._opright)
            
        if args:
            args = "{},{}".format(var,args)
        else:
            args = var
                
        return "{}({})".format(self._op_name,args)
        
    def __repr__(self):
        if self._opright is not None:
            opright = " * {}".format(self._opright.__repr__())
        else:
            opright = ""
            
        return "{}({}){}".format(self.__class__.__name__,self._op_args,opright)
        
    def __mul__(self,rother):
        # Op1 * Op2
        if rother is None or isinstance(rother,Identity):
            return self
        elif isinstance(rother,Operator):
            op = copy.copy(self)
            if self._opright is None:
                op._opright = rother
            else:
                op._opright = self._opright * rother
            return op
        else:
            raise NotImplementedError

    def __rmul__(self,lother):
        # first self, then lother
        if isinstance(lother,Operator): 
            return lother.__mul__(self)
        elif lother is None:
            return self
        else:
            raise NotImplementedError
            
    def __call__(self,x):
        if self._opright is not None:
            x = self._opright(x)
        return self._eval(x)

    def _eval(self,arg):
        raise NotImplementedError
        
    @property
    def inverse(self):
        if self._opright is None:
            return self._inverse
        else:
            return self._opright.inverse * self._inverse

    @property
    def _inverse(self):
        raise NotImplementedError

    def __eq__(self,other):
        if not isinstance(other,self.__class__):
            return False
        return self._op_name==other._op_name and self._attr_eq(other)

    def _attr_eq(self,other):
        return True    
    
    def __ne__(self,other):
        return not self.__eq__(other)

    def tofloat(self,x):
        # Because x could be a pint.Quantity
        return x*np.float64(1)


class Identity(Operator):

    @property
    def _op_name(self):
        return "id"

    def __mul__(self,rother):
        if rother is None:
            return self  
        else:
            return rother
        
    def _eval(self,arg):
        return arg

    @property
    def _inverse(self):
        return self


class Lambda(Operator):

    def __init__(self,name,iname,func,ifunc):
        self.name,self.iname,self.func,self.ifunc = name,iname,func,ifunc
        super(Lambda,self).__init__()

    @property
    def _op_name(self):
        return self.name
        
    def _eval(self,arg):
        return self.func(arg)
        
    @property
    def _inverse(self):
        return self.__class__(self.iname,self.name,self.ifunc,self.func)


class CommutativeOperator(Operator):
    """ f: x -> f(x): f*LinearOperator = LinearOperator*g and g exists
    """
    
    def __mul__(self,rother):
        if isinstance(rother,LinearOperator):
            if self._opright is not None:
                return super(CommutativeOperator,self).__mul__(rother)
            op = copy.copy(rother)
            op._opright = self._commute_linop(rother) * rother._opright
            return op   
        else:
            return super(CommutativeOperator,self).__mul__(rother)

    def _commute_linop(self,linop):
        """If operator f is self, get operator g for which f * linop = linop * g
        
        Args:
            linop(LinearOperator):
            
        Returns:
            CommutativeOperator: g
        """
        raise NotImplementedError
        

class ClipOperator(CommutativeOperator):
    """ f: x -> x       cmin < x < cmax
                ...     else
    """
    
    def __init__(self,cmin=None,cmax=None):
        self.cmin = cmin
        self.cmax = cmax
        super(ClipOperator,self).__init__()

    def _attr_eq(self,other):
        return self.cmin==other.cmin and self.cmax==other.cmax
        
    @property
    def _op_args(self):
        return "{},{}".format(self.cmin,self.cmax)

    def _commute_linop(self,linop):
        m = self.tofloat(linop.m)
        with np.errstate(divide='ignore', invalid='ignore'):
            if self._valid_limit(self.cmin):
                cmin = (self.cmin-linop.b)/m
            else:
                cmin = None
            if self._valid_limit(self.cmax):
                cmax = (self.cmax-linop.b)/m
            else:
                cmax = None
            if units.binary_operator(m,0,operator.lt):
                cmin,cmax = cmax,cmin
        return self.__class__(cmin,cmax)

    def _valid_limit(self,c):
        if c is None:
            return False
        return not np.isnan(c)
        
    def _valid_limits(self,cmin,cmax):
        vmin = self._valid_limit(cmin)
        vmax = self._valid_limit(cmax)
        
        if vmin and vmax:
            v = self.cmin<=self.cmax
        else:
            v = True
        
        return v,vmin,vmax
        
    def __mul__(self,rother):
        if isinstance(rother,self.__class__):
            if self._opright is not None:
                return super(ClipOperator,self).__mul__(rother)

            vs,vsmin,vsmax = self._valid_limits(self.cmin,self.cmax)
            vr,vrmin,vrmax = self._valid_limits(rother.cmin,rother.cmax)
            
            bnocomb = not vs or not vr
            bnocomb |= vsmax and vrmin and self.cmax<rother.cmin
            bnocomb |= vrmax and vsmin and rother.cmax<self.cmin
            if bnocomb:
                return super(ClipOperator,self).__mul__(rother)

            cmin,cmax = self.cmin,self.cmax
            if vrmin:
                if vsmin:
                    cmin = max(cmin,rother.cmin)
                else:
                    cmin = rother.cmin
                    
            if vrmax:
                if vsmax:
                    cmax = min(cmax,rother.cmax)
                else:
                    cmax = rother.cmax

            op = self.__class__(cmin,cmax)
            op._opright = rother._opright
            return op
        else:
            return super(ClipOperator,self).__mul__(rother)


class Clip(ClipOperator):
    """ f: x -> x       cmin < x < cmax
                cmin    cmin < x
                cmax    x < cmax
    """
    
    # Expressed with Heaviside Step Function:
    # min(x,cmax) = cmax.H(x-cmax) + x.[1-H(x-cmax)]
    # max(x,cmin) = cmin.[1-H(x-cmin)] + x.H(x-cmin)
    # clip(x,cmin,cmax) = min(max(x,cmin),cmax)
    
    @property
    def _op_name(self):
        return "clip"
    
    def _eval(self,x):
        #return np.clip(x,self.cmin,self.cmax)
        y,func = units.asqarrayf(x)
        y = self.tofloat(y)
        
        v,vmin,vmax = self._valid_limits(self.cmin,self.cmax)
        if v:
            if vmin:
                y[units.binary_operator(y,self.cmin,operator.lt)] = self.cmin
            if vmax:
                y[units.binary_operator(y,self.cmax,operator.gt)] = self.cmax
        else:
            nan = units.quantity_like(np.nan,y)
            y[:] = nan

        return func(y)

    @property
    def _inverse(self):
        return NaNClip(self.cmin,self.cmax)

        
class NaNClip(ClipOperator):
    """ f: x -> x    cmin < x < cmax
                nan  else
    """
    
    @property
    def _op_name(self):
        return "nclip"

    def _eval(self,x):
        y,func = units.asqarrayf(x)
        y = self.tofloat(y)

        nan = units.quantity_like(np.nan,y)
        if self._valid_limit(self.cmin):
            y[units.binary_operator(y,self.cmin,operator.lt)] = nan
        if self._valid_limit(self.cmax):
            y[units.binary_operator(y,self.cmax,operator.gt)] = nan
        
        return func(y)

    @property
    def _inverse(self):
        return NaNClip(self.cmin,self.cmax)


class LinearOperator(Operator):
    """ f: x -> m*x+b
    """
    def __init__(self,m,b):
        self.m = m
        self.b = b
        super(LinearOperator,self).__init__()
        
    @property
    def _op_args(self):
        return "{},{}".format(self.m,self.b)
        
    def __str__(self):
        if self._opright is None:
            x = "x"
        else:
            x = str(self._opright)
        mx = "{} * {}".format(self.m,x)

        b = " {:+}".format(self.b)
            
        return "{}{}".format(mx,b)
        
    def _attr_eq(self,other):
        return self.m==other.m and self.b==other.b
        
    def _eval(self,x):
        return self.m*x + self.b

    def __mul__(self,rother):
        if isinstance(rother,self.__class__):
            if self._opright is not None:
                return super(LinearOperator,self).__mul__(rother)
                
            op = self.__class__(self.m*rother.m,rother.b*self.m+self.b)
            op._opright = rother._opright
            return op
        else:
            return super(LinearOperator,self).__mul__(rother)

    def __pow__(self,p):
        if isinstance(p,numbers.Integral):
            if p<0:
                raise NotImplementedError
            elif p==0:
                ret = self.__class__(0,1)
            else:
                ret = self.__class__(self.m,self.b)
                for i in range(1,p):
                    ret *= self.__class__(self.m,self.b)
            
            return ret
        else:
            raise NotImplementedError

    @property
    def inverse(self):
        m = self.tofloat(self.m)
        with np.errstate(divide='ignore', invalid='ignore'):
            selfi = self.__class__(1/m,-self.b/m)
        if self._opright is None:
            return selfi
        else:
            return self._opright.inverse * selfi




