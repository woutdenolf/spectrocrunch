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

import pint
ureg = pint.UnitRegistry()
Q_ = ureg.Quantity

def speedoflight():
    return Q_(299792458, ureg.m/ureg.s)

def planckconstant():
    return Q_(4.13566743E-15, ureg.eV*ureg.s)

def wavelengthenergy(x,keV=False):
    # E(keV) = h(eV sec) . c(m/sec)
    #          ------------------
    #           lambda(nm).10^-9
    
    if not isinstance(x,Q_):
        x = Q_(x, ureg.nm)
    x.ito(ureg.m)
    ret = (planckconstant()*speedoflight())/x
    if kev:
        ret.ito(ureg.keV)
    return ret

def elementarycharge():
    return Q_(1.60217662E-19, ureg.C) # C or J/eV

def temperatureinkelvin(T):
    if not isinstance(T,Q_):
        T = Q_(T,ureg.degC)
    return T.to(ureg.kelvin)
       
def eholepair_si(T=21):
    # https://doi.org/10.1016/j.nima.2007.03.020
    T = temperatureinkelvin(T)
    x = [80,270] * ureg.kelvin
    y = [3.77,3.68] * ureg.eV
    
    m = (y[1]-y[0])/(x[1]-x[0])
    b = y[1]-m*x[1]
    
    return m*T+b # eV
    
