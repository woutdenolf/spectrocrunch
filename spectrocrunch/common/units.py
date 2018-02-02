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
        return x
    except:
        return ureg.Quantity(x,units=units)

def units(x):
    try:
        return x.units
    except:
        return None

def generator(x):
    kwargs = {"magnitude":x.magnitude,"units":x.units}
    return {"generator":QuantityGenerator,"kwargs":kwargs}

def magnitude(x,units=None):
    """Magnitude when Quantity, untouched otherwise.
    """
    try:
        return x.to(units).magnitude
    except:
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
    except:
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
    except:
        return x
        


