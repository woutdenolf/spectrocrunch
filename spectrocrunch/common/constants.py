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
from ..common import units

def temperatureinkelvin(T):
    """
    Args:
        T(num|array): temperature in deg C
    Returns:
        num|array: keV
    """
    T = units.Quantity(T,ureg.degC)
    return T.to(ureg.kelvin)
       
def eholepair_si(T=21):
    """
    Args:
        T(num|array): temperature in deg C
    Returns:
        num|array: keV
    """
    # https://doi.org/10.1016/j.nima.2007.03.020
    T = temperatureinkelvin(T)
    x = units.Quantity([80,270],ureg.kelvin) # K
    y = units.Quantity([3.77,3.68],"eV")
    
    m = (y[1]-y[0])/(x[1]-x[0])
    b = y[1]-m*x[1]
    
    return (m*T+b).to("eV")
    
    