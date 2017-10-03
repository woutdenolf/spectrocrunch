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

import unittest

from ..compound import compound as compoundraw
from ..compoundfromformula import compoundfromformula
from ..compoundfromlist import compoundfromlist
from ..compoundfromcif import compoundfromcif
from ..compoundfromname import compoundfromname as compoundfromname
from ..types import fractionType
from ..element import element
from ... import ureg

try:
    import iotbx.cif as iotbxcif
except ImportError:
    iotbxcif = None

import numpy as np
import xraylib

class test_compound(unittest.TestCase):

    def test_comparable(self):
        c1 = compoundfromformula("C6H2(NO2)3CH3",1.2,name="compound")
        c2 = compoundfromformula("C6H2(NO2)3CH3",1.2,name="compound")
        c3 = compoundfromformula("C6H2(NO2)3CH3",1.21,name="compound")
        c4 = compoundfromformula("C6H2(NO2)3CH3",1.2,name="compoundx")
        c5 = compoundfromformula("C6H2(NO2)3CH2",1.2,name="compound")
        self.assertEqual(c1,c2)
        self.assertEqual(c1,c3)
        self.assertNotEqual(c1,c4)
        self.assertEqual(c1,c5) # this is by design but may be unwanted?

        self.assertEqual(element("Ca"),element("Ca"))
        self.assertNotEqual(element("Ca"),element("C"))
        self.assertEqual(element("Ca"),"Ca")
        self.assertNotEqual(element("Ca"),"C")

    def test_formula(self):
        elements = ["C","N","O","H"]
        a = [7,3,6,5]
        density = 2.3
        c = compoundfromformula("C6H2(NO2)3CH3",density)

        elements2 = c.molefractions()
        for i in range(len(elements)):
            self.assertEqual(elements2[elements[i]],a[i])

        self.assertEqual(c.density,density)

    def test_list(self):
        elements = ["Fe","S","O"]
        a = [1,1,4.]
        density = 2.3
        c = compoundfromlist(elements,a,fractionType.mole,density)

        elements2 = c.molefractions()
        for i in range(len(elements)):
            self.assertEqual(elements2[elements[i]],a[i])

        self.assertEqual(c.density,density)

        c = compoundfromlist(elements,a,fractionType.weight,density)
        wfrac = c.weightfractions()
        for i in range(len(elements)):
            self.assertAlmostEqual(wfrac[elements[i]],a[i]/float(sum(a)))

    def test_cif(self):
        if iotbxcif is None:
            raise unittest.SkipTest("cctbx not available")

        elements = ["Ca","C","O"]
        a = [6,6,18.] # unit cell content
        c = compoundfromcif("cif/calcite.cif",name="calcite") 

        elements2 = c.molefractions()
        for i in range(len(elements)):
            self.assertEqual(elements2[elements[i]],a[i])

    def test_addelement(self):
        c1 = compoundraw(["Fe","O"],[2,3],fractionType.mole,density=1)
        c2 = compoundraw(["Fe"],[2],fractionType.mole,density=1)
        c2.addelement("O",3,fractionType.mole)
        self.assertEqual(c1.elements,c2.elements)

        c1 = compoundraw(["Fe","O"],[2,3],fractionType.weight,density=1)
        c2 = compoundraw(["Fe"],[2],fractionType.weight,density=1)
        c2.addelement("O",3/5.,fractionType.weight)
        self.assertEqual(c1.elements.keys(),c2.elements.keys())
        for v1,v2 in zip(c1.elements.values(),c2.elements.values()):
            self.assertAlmostEqual(v1,v2)

    def test_vacuum(self):
        c = compoundraw([],[],fractionType.mole,0,name="vacuum")
        self.assertEqual(len(c.elements),0)
        self.assertEqual(len(c.weightfractions()),0)
        self.assertEqual(len(c.molefractions()),0)
        self.assertEqual(c.molarmass(),0)
        self.assertEqual(c.density,0)

    def test_name(self):
        c = compoundfromname("vacuum")
        self.assertEqual(len(c.elements),0)
        self.assertEqual(len(c.weightfractions()),0)
        self.assertEqual(len(c.molefractions()),0)
        self.assertEqual(c.molarmass(),0)
        self.assertEqual(c.density,0)

        c = compoundfromname("linseed oil")

        with self.assertRaises(KeyError) as context:
            c = compoundfromname("linseed oill")

    def test_refractiveindex(self):
        density = 2.328
        c = compoundfromformula("Si",density,name="calcite")
        energy = np.asarray([5,10])
        
        Z = 14
        
        # Xraylib
        n_re0 = np.asarray([xraylib.Refractive_Index_Re("Si",e,density) for e in energy])
        n_im0 = np.asarray([xraylib.Refractive_Index_Im("Si",e,density) for e in energy])
        
        # Exactly like Xraylib
        n_re1 = 1-density*4.15179082788e-4*(Z+np.asarray([xraylib.Fi(Z, e) for e in energy]))/xraylib.AtomicWeight(Z)/energy**2;
        n_im1 = np.asarray([xraylib.CS_Total(Z, e) for e in energy])*density*9.8663479e-9/energy
        
        np.testing.assert_allclose(n_re0,n_re1)
        np.testing.assert_allclose(n_im0,n_im1)
        
        # Im: Kissel
        n_im1b = np.asarray([xraylib.CS_Total_Kissel(Z, e) for e in energy])*density*9.8663479e-9/energy
        
        # Im: Kissel with pint
        n_im1d = ureg.Quantity(c.mass_att_coeff(energy),'cm^2/g') *\
                ureg.Quantity(c.density,'g/cm^3') *\
                ureg.Quantity(energy,'keV').to("cm",'spectroscopy')/(4*np.pi)
        n_im1d = n_im1d.to('dimensionless').magnitude
        
        r = n_im1b/n_im1d
        np.testing.assert_allclose(r,r[0])

        # Im: other formula
        n_im1c = -density*4.15179082788e-4*(np.asarray([xraylib.Fii(Z, e) for e in energy]))/xraylib.AtomicWeight(Z)/energy**2;

        # Spectrocrunch
        n_re2 = c.refractive_index_re(energy)
        n_im2 = c.refractive_index_im(energy)

        np.testing.assert_allclose(n_re2,n_re0)
        np.testing.assert_allclose(n_im2,n_im1c,rtol=1e-6)
        
def test_suite_all():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_compound("test_comparable"))
    testSuite.addTest(test_compound("test_formula"))
    testSuite.addTest(test_compound("test_list"))
    testSuite.addTest(test_compound("test_cif"))
    testSuite.addTest(test_compound("test_addelement"))
    testSuite.addTest(test_compound("test_vacuum"))
    testSuite.addTest(test_compound("test_name"))
    testSuite.addTest(test_compound("test_refractiveindex"))
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite_all()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
