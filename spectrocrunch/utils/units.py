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

from ..patch.pint import ureg
from ..patch.pint import errors as pinterrors
from . import instance
from . import listtools

import numpy as np
from jsonpickle.handlers import BaseHandler,register

def Quantity(x,units=None,forcequantity=True):
    """Quantity with given units when not present
    
    Args:
        x(num|array|Quantity):
        units(Optional(Unit)): dimensionless by default
        forcequantity(Optional(bool)): returns Quantity when forcequantity=True
                
    returns:
        x(num|array|Quantity): magnitude is scalar or array, depending on the input
    """
    if instance.isarray(x):
        x = asqarray(x,forcequantity=False)
    if forcequantity:
        if not instance.isquantity(x):
            x = ureg.Quantity(x,units=units)
    return x
    
def magnitude(x,units=None):
    """Magnitude of a Quantity or number (units not ignored when number)
    
    Args:
        x(num|array|Quantity):
        units(Optional(Unit)): dimensionless by default
    
    Returns:
        x(num|array)
    
    Raises:
        DimensionalityError
    """
    if units is None:
        units = ureg.dimensionless
    
    if instance.isquantity(x):
        return x.to(units).magnitude
        
    if instance.isarraynot0(x):
        return np.asarray([magnitude(y,units=units) for y in x])
        #TODO: not sure why this fails:
        #return np.vectorize(lambda y: magnitude(y,units=units))(x) 
    else:
        if units != ureg.dimensionless:
            raise pinterrors.DimensionalityError(None,units)
        return x

def umagnitude(x,units=None):
    """Magnitude of a Quantity or number (units ignored when number)
    
    Args:
        x(num|array|Quantity):
        units(Optional(Unit)): dimensionless by default
    
    Returns:
        x(num|array)
    """
    return Quantity(x,units).to(units).magnitude

def units(x):
    """Units of a Quantity
    
    Args:
        x(num|array|Quantity): 
        
    Returns:
        x(None|Unit|array(Unit)): 
    """
    
    if instance.isquantity(x):
        return x.units
        
    if instance.isarraynot0(x):
        return [units(y) for y in x]
    else:
        return None

def asqarray(x,forcequantity=True,**kwargs):
    """Quantity with array magnitude
    
    Args:
        x(num|array|Quantity):
        forcequantity(Optional(bool)): returns Quantity when forcequantity=True
    
    Returns:
        x(Quantity|array):
    """
    
    u = next(listtools.flatten(units(x)))
    
    if instance.isarraynot0(x):
        x = [magnitude(y,units=u) for y in x]
    else:
        x = magnitude(x,u)

    x = instance.asarray(x,**kwargs)
    
    if forcequantity or u is not None:
        x = Quantity(x,units=u)

    return x
    
def asqarrayf(x,**kwargs):
    """Quantity with array magnitude, including scalar restoration
        
    Args:
        x(num|array|Quantity):

    Returns:
        x(Quantity|array):
        func(callable): apply to x to restore scalars
    """
    if instance.isarray(x):
        func = lambda x:x[0]
    else:
        func = lambda x:x

    return asqarray(x,**kwargs),func

def quantity_like(x,y,forcequantity=True):
    """Quantity with units of y
        
    Args:
        x(num|array|Quantity):
        y(num|array|Quantity):
        forcequantity(Optional(bool)): returns Quantity when forcequantity=True
    
    Returns:
        x(Quantity):
    """
    u = set(listtools.flatten(units(y)))
    if len(u)!=1:
        raise ValueError("Array values have different units")
    u = next(iter(u))
    if not forcequantity and u is None:
        return magnitude(x,units=None)
    else:
        if u is None:
            u = ureg.dimensionless
        return Quantity(x,units=u).to(u)
    
def unitsto(x,units):
    """Converts to the first valid unit
    
    Args:
        x(array(Quantity)|Quantity):
        units(Unit|array(Unit)):
    
    Returns:
        x(Quantity)
    """
    if not instance.isarray(units):
        units = [units]
    x = Quantity(x)
    for unit in units:
        try:
            x = x.to(unit)
            break
        except pinterrors.DimensionalityError:
            continue
    else:
        raise RuntimeError("Units of {} cannot be converted to {}".format(x,units))
    return x

def binary_operator(a,b,op):
    """Compare quantities and numbers
    
    Args:
        x(num|Quantity):
        y(num|Quantity):
        op(callable):

    Returns:
        op(a,b)
    """
    return op(Quantity(a,forcequantity=False),quantity_like(b,a,forcequantity=False))

class QuantityHandler(BaseHandler):
    
    def flatten(self, quantity, data):
        data['magnitude'] = quantity.magnitude
        data['units'] = str(quantity.units)
        return data
    
    def restore(self, data):
        return ureg.Quantity(data['magnitude'],units=data['units'])

def jsonpickle_register_handlers():
    register(ureg.Quantity, QuantityHandler, base=True)
