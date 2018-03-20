# -*- coding: utf-8 -*-
#
#   Copyright (C) 2015 European Synchrotron Radiation Facility, Grenoble, France
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

execfile("initcctbx.py")


# Don't use the installed version
import os, sys
sys.path.insert(1,os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from spectrocrunch.materials.compoundfromformula import compoundfromformula

from spectrocrunch.materials.mixture import mixture
from spectrocrunch.materials.types import fraction

import scipy.integrate as integrate
import scipy.optimize as optimize

import numpy as np

def test1():
    compound1 = compoundfromformula("La2O3",0,name="La2O3")
    compound2 = compoundfromformula("SrO",0,name="SrO")
    compound3 = compoundfromformula("Co2O3",0,name="Co2O3")
    compound4 = compoundfromformula("Fe2O3",0,name="Fe2O3")

    m = mixture([compound1,compound2,compound3,compound4],[1,1,1,1],fraction.mass)

    print(compound4.massfractions())
    print("")
    elements = m.elemental_massfractions()
    print("")
    for e in elements:
        print(e,e.MM,elements[e])

class capillary_transmission:

    def __init__(self,mu,rho,packing=1.):
        self.mu = mu
        self.rho = rho*packing

    def __call__(self,R):
        if R == 0:
            return 1.
        else:
            #return np.exp(-2.*R*self.mu*self.rho) # estimation
            return integrate.quad(lambda x: np.exp(-2*self.mu*self.rho*np.sqrt(R*R-x*x)), -R, R)[0]/(2*R)

class capillary_transmission2:

    def __init__(self,mu,rho,R):
        self.mu = mu
        self.rho = rho
        self.R = R

    def __call__(self,packing):
        return integrate.quad(lambda x: np.exp(-2*self.mu*self.rho*packing*np.sqrt(self.R*self.R-x*x)), -self.R, self.R)[0]/(2*self.R)

class capillary_refine:

    def __init__(self,mu,rho,T,packing):
        self.o = capillary_transmission(mu,rho,packing=packing)
        self.transmission = T

    def __call__(self,R):
        return self.o(R)-self.transmission

class capillary_refine2:

    def __init__(self,mu,rho,T,R):
        self.o = capillary_transmission2(mu,rho,R)
        self.transmission = T

    def __call__(self,R):
        return self.o(R)-self.transmission


def test2():
    compound1 = compoundfromformula("C18H30O2",0.9291,name="linseedoil")
    compound2 = compoundfromformula("Pb3C2O8H2",6.8,name="hydrocerussite")
    m = mixture([compound1,compound2],[0.5,0.5],fraction.volume)
    mu = m.mass_att_coeff(35.)

    # 60% transmission
    T = 0.6 

    # flat sample
    thickness = -np.log(T)/(mu*m.density) # cm

    # capillary
    #R = 30e-4 # cm
    #packing = optimize.broyden2(capillary_refine2(mu,m.density,T,R),0.7)
    packing = 1.
    R = optimize.broyden2(capillary_refine(mu,m.density,T,packing),[80e-4])[0] # cm

    # 1-3 ug
    mass = 1*1e-6 # g
    volume = mass/m.density # cm^3

    print "Mixture:"
    print "density = {} g/cm^3".format(m.density)
    print "mass.att. = {} cm^2/g".format(mu)

    print "\nCapillary:"
    print "R = {} um".format(R*1e4)
    print "packing = {} %".format(packing*100)
    print "h = {} mm (@ total mass = {} ug)".format(volume/(np.pi*R*R)*10,mass*1e6)


    print "\nFlat sample:"
    print "thickness = {} um".format(thickness*1e4)
    print "footprint = {} mm^2 (@ total mass = {} ug)".format(volume/thickness*1e4,mass*1e6)
    

if __name__ == '__main__':
    test2()


    

    

