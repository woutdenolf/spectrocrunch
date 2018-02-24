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

from .. import compound as compoundraw
from .. import compoundfromformula
from .. import compoundfromlist
from .. import compoundfromcif
from .. import compoundfromname
from ..types import fractionType
from ... import ureg

try:
    import iotbx.cif as iotbxcif
except ImportError:
    iotbxcif = None

import numpy as np
import xraylib

class test_compound(unittest.TestCase):

    def test_comparable(self):
        c1 = compoundfromformula.CompoundFromFormula("C6H2(NO2)3CH3",1.2,name="compound")
        c2 = compoundfromformula.CompoundFromFormula("C6H2(NO2)3CH3",1.2,name="compound")
        c3 = compoundfromformula.CompoundFromFormula("C6H2(NO2)3CH3",1.21,name="compound")
        c4 = compoundfromformula.CompoundFromFormula("C6H2(NO2)3CH3",1.2,name="compoundx")
        c5 = compoundfromformula.CompoundFromFormula("C6H2(NO2)3CH2",1.2,name="compound")
        self.assertEqual(c1,c2)
        self.assertEqual(c1,c3)
        self.assertNotEqual(c1,c4)
        self.assertEqual(c1,c5) # this is by design but may be unwanted?

    def test_formula(self):
        elements = ["C","N","O","H"]
        a = [7,3,6,5]
        density = 2.3
        c = compoundfromformula.CompoundFromFormula("C6H2(NO2)3CH3",density)

        elements2 = c.molefractions()
        for i in range(len(elements)):
            self.assertEqual(elements2[elements[i]],a[i])

        self.assertEqual(c.density,density)

    def test_list(self):
        elements = ["Fe","S","O"]
        a = [1,1,4.]
        density = 2.3
        c = compoundfromlist.CompoundFromList(elements,a,fractionType.mole,density)

        elements2 = c.molefractions()
        for i in range(len(elements)):
            self.assertEqual(elements2[elements[i]],a[i])

        self.assertEqual(c.density,density)

        c = compoundfromlist.CompoundFromList(elements,a,fractionType.weight,density)
        wfrac = c.weightfractions()
        for i in range(len(elements)):
            self.assertAlmostEqual(wfrac[elements[i]],a[i]/float(sum(a)))

    def test_cif(self):
        if iotbxcif is None:
            raise unittest.SkipTest("cctbx not available")

        elements = ["Ca","C","O"]
        a = [6,6,18.] # unit cell content
        c = compoundfromcif.CompoundFromCif("cif/calcite.cif",name="calcite") 

        elements2 = c.molefractions()
        for i in range(len(elements)):
            self.assertEqual(elements2[elements[i]],a[i])

    def test_addelement(self):
        c1 = compoundraw.Compound(["Fe","O"],[2,3],fractionType.mole,density=1)
        c2 = compoundraw.Compound(["Fe"],[2],fractionType.mole,density=1)
        c2.addelement("O",3,fractionType.mole)
        self.assertEqual(c1.elements,c2.elements)

        c1 = compoundraw.Compound(["Fe","O"],[2,3],fractionType.weight,density=1)
        c2 = compoundraw.Compound(["Fe"],[2],fractionType.weight,density=1)
        c2.addelement("O",3/5.,fractionType.weight)
        n1 = c1.molefractions()
        n2 = c2.molefractions()
        self.assertEqual(sorted(n1.keys()),sorted(n2.keys()))
        for k in n1:
            self.assertAlmostEqual(n1[k],n2[k])

    def test_vacuum(self):
        c = compoundraw.Compound([],[],fractionType.mole,0,name="vacuum")
        self.assertEqual(len(c.elements),0)
        self.assertEqual(len(c.weightfractions()),0)
        self.assertEqual(len(c.molefractions()),0)
        self.assertEqual(c.molarmass(),0)
        self.assertEqual(c.density,0)

    def test_name(self):
        c = compoundfromname.compoundfromname("vacuum")
        self.assertEqual(len(c.elements),0)
        self.assertEqual(len(c.weightfractions()),0)
        self.assertEqual(len(c.molefractions()),0)
        self.assertEqual(c.molarmass(),0)
        self.assertEqual(c.density,0)

        c = compoundfromname.compoundfromname("linseed oil")

        with self.assertRaises(KeyError) as context:
            c = compoundfromname.compoundfromname("linseed oill")

    def test_refractiveindex(self):
        density = 2.328
        c = compoundfromformula.CompoundFromFormula("Si",density,name="silicon")
        energy = np.asarray([5,10])
        
        Z = 14
        
        # Xraylib (TODO: n_im sign is wrong!)
        n_re0 = np.asarray([xraylib.Refractive_Index_Re("Si",e,density) for e in energy])
        n_im0 = -np.asarray([xraylib.Refractive_Index_Im("Si",e,density) for e in energy])
        
        # Exactly like Xraylib (TODO: CS_Total -> CS_Photo_Total)
        delta = density*4.15179082788e-4*(Z+np.asarray([xraylib.Fi(Z, e) for e in energy]))/xraylib.AtomicWeight(Z)/energy**2
        beta = np.asarray([xraylib.CS_Total(Z, e) for e in energy])*density*9.8663479e-9/energy
        n_re1 = 1-delta
        n_im1 = -beta
        
        np.testing.assert_allclose(n_re0,n_re1)
        np.testing.assert_allclose(n_im0,n_im1)
        
        # Im: Kissel
        n_im1b = -np.asarray([xraylib.CS_Total_Kissel(Z, e) for e in energy])*density*9.8663479e-9/energy
        
        # Im: Kissel with pint
        n_im1d = ureg.Quantity(c.mass_att_coeff(energy),'cm^2/g') *\
                ureg.Quantity(c.density,'g/cm^3') *\
                ureg.Quantity(energy,'keV').to("cm",'spectroscopy')/(4*np.pi)
        n_im1d = -n_im1d.to('dimensionless').magnitude
        
        r = n_im1b/n_im1d
        np.testing.assert_allclose(r,r[0])

        # Im: other formula
        n_im1c = density*4.15179082788e-4*(np.asarray([xraylib.Fii(Z, e) for e in energy]))/xraylib.AtomicWeight(Z)/energy**2;

        # Spectrocrunch
        n_re2 = c.refractive_index_re(energy)
        n_im2 = c.refractive_index_im(energy)

        np.testing.assert_allclose(n_re2,n_re0)
        np.testing.assert_allclose(n_im2,n_im1c,rtol=1e-6)
    
    def test_refractiveindex2(self):
        density = 5.3
        c = compoundfromformula.CompoundFromFormula("Fe2O3",5.3,name="test")
        energy = 30
        
        wavelength = ureg.Quantity(energy,'keV').to("cm","spectroscopy")
        beta = c.refractive_index_beta(energy)
        m = (4*np.pi/(wavelength*ureg.Quantity(c.density,'g/cm^3'))).to("cm^2/g").magnitude
        np.testing.assert_allclose(beta*m,c.mass_abs_coeff(energy),rtol=1e-2)
        
        delta = c.refractive_index_delta(energy)
        m = (2*np.pi/(ureg.re*wavelength**2*ureg.avogadro_number/ureg.Quantity(c.molarmasseff(),"g/mol")*c.Zeff)).to("g/cm^3").magnitude
        np.testing.assert_allclose(delta*m,c.density,rtol=1e-2)
        
        
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
    testSuite.addTest(test_compound("test_refractiveindex2"))
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite_all()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
