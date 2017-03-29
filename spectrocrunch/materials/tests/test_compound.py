# -*- coding: utf-8 -*-
#
#   Copyright (C) 2016 European Synchrotron Radiation Facility, Grenoble, France
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

from ..compoundfromformula import compoundfromformula
from ..compoundfromlist import compoundfromlist
from ..compoundfromcif import compoundfromcif
from ..types import fractionType
from ..element import element

try:
    import iotbx.cif as iotbxcif
except ImportError:
    iotbxcif = None

import numpy as np

class test_compound(unittest.TestCase):

    def test_comparable(self):
        c1 = compoundfromformula("C6H2(NO2)3CH3",1.2,name="compound")
        c2 = compoundfromformula("C6H2(NO2)3CH3",1.2,name="compound")
        c3 = compoundfromformula("C6H2(NO2)3CH3",1.21,name="compound")
        c4 = compoundfromformula("C6H2(NO2)3CH3",1.2,name="compoundx")
        self.assertEqual(c1,c2)
        self.assertEqual(c1,c3)
        self.assertNotEqual(c1,c4)

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
        a = [6,6,18.]
        c = compoundfromcif("calcite",name="calcite") 

        elements2 = c.molefractions()
        for i in range(len(elements)):
            self.assertEqual(elements2[elements[i]],a[i])

def test_suite_all():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_compound("test_comparable"))
    testSuite.addTest(test_compound("test_formula"))
    testSuite.addTest(test_compound("test_list"))
    testSuite.addTest(test_compound("test_cif"))
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite_all()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
