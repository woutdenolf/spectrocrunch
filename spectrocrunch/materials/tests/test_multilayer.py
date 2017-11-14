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

from .. import multilayer
from .. import compoundfromformula
from .. import mixture
from .. import types
from .. import utils

import numpy as np

class test_compound(unittest.TestCase):

    def _multilayer(self):
        c1 = compoundfromformula.CompoundFromFormula("CaCO3",2.71,name="calcite")
        c2 = compoundfromformula.CompoundFromFormula("Fe2O3",5.3,name="hematite")
        c3 = compoundfromformula.CompoundFromFormula("PbCO3",6.53,name="cerussite")

        l = [c1,mixture.Mixture([c2,c3],[1,1],types.fractionType.mole)]
        thickness = [10,20]
        o = multilayer.Multilayer(material=l, thickness=thickness, anglein=45, angleout=45+90, solidangle=4*np.pi*0.1)
        
        return o,thickness
        
    def test_str(self):
        o,thickness = self._multilayer()
        s = str(o).split("\n")
        for i,layer in enumerate(o):
            self.assertTrue(str(layer) in s[i+1])
            self.assertTrue(str(thickness[i]) in s[i+1])
            
    def test_transmission(self):
        o,thickness = self._multilayer()
        
        energy = np.linspace(7,8,5)

        T1 = o.transmission(energy)
        T2 = o.transmission(energy,decomposed=True)
        
        T3 = np.ones_like(energy)
        for mat in T2:
            for cs in utils.elemental_crosssections(mat["cs"]).values():
                T3 *= np.exp(-mat["density"]*mat["thickness"]*cs)

        np.testing.assert_allclose(T1,T3)
        
def test_suite_all():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_compound("test_str"))
    testSuite.addTest(test_compound("test_transmission"))
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite_all()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
