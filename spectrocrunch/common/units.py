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

from .. import ureg
from . import instance
from . import persistence

def Quantity(x,units=None):
    """Add units when not present. Use instead of ureg.Quantity in case x may already be a Quantity.
    """
    try:
        x.units
    except AttributeError:
        return ureg.Quantity(x,units=units)
    else:
        return x
        
def units(x):
    try:
        return x.units
    except AttributeError:
        return None

def to(x,y):
    """Convert units of x to units of y (if y has no units, it will try to have x without units as we)
    """
    try:
        y.units
    except AttributeError: # y has no units
        try:
            x.units 
        except AttributeError: # x has no units
            return x
        else:
            try:
                x = x.to("dimensionless").magnitude
            except:
                raise RuntimeError("Units are not compatible: {}, {}".format(x,y))
    else: # y has units
        try:
            x = x.to(y.units)
        except:
            raise RuntimeError("Units are not compatible: {}, {}".format(x,y))
            
    return x

def generator(x):
    kwargs = {"magnitude":x.magnitude,"units":x.units}
    return {"generator":QuantityGenerator,"kwargs":kwargs}

def magnitude(x,units=None):
    """Magnitude when Quantity, untouched otherwise.
    """
    try:
        return x.to(units).magnitude
    except AttributeError:
        return x

def quantity_like(x,y):
    try:
        return Quantity(x,units=y.units)
    except:
        return x

def binary_operator(a,b,op):
    """Compare quantities and numbers
    """
    try:
        return op(a,b)
    except:
        return op(magnitude(a),magnitude(a))

def asarrayf(x,**kwargs):
    u = units(x)
    if u is None:
        funcu = lambda x:x
    else:
        funcu = lambda x:Quantity(x,units=u)
        
    try:
        x = x.magnitude
    except AttributeError:
        pass

    x,func = instance.asarrayf(x,**kwargs)

    funcw = lambda x: funcu(func(x))

    return funcu(x),funcw

def serialize(x):
    try:
        x.units
        return persistence.SerializedGenerator(module=__name__,\
                    generator="Quantity",\
                    args=(x.magnitude,),\
                    kwargs={"units":str(x.units)})
    except AttributeError:
        return x
        


