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

from ..linop import linop

class test_ops(unittest.TestCase):

    def test_linop(self):
        for x in range(-2,2):
            for m1 in range(-11,11,3):
                for b1 in range(-11,11,3):
                    o1 = linop(m1,b1)
                    o1i = o1**(-1)

                    self.assertAlmostEqual(o1(x),o1.m*x+o1.b)
                    self.assertAlmostEqual(o1i(x),(x-o1.b)/float(o1.m))
                    self.assertEqual(o1*o1,o1**2)
                    self.assertEqual(o1*o1*o1,o1**3)
                    self.assertAlmostEqual((o1i*o1i)(x),(o1**(-2))(x))
                    self.assertAlmostEqual((o1i*o1i*o1i)(x),(o1**(-3))(x))

                    for m2 in range(-11,11,3):
                        for b2 in range(-11,11,3):
                            o2 = linop(m2,b2)
                            o2i = o2**(-1)
                            self.assertAlmostEqual((o1*o2)(x),o2(o1(x)))
                            self.assertAlmostEqual((o1/o2)(x),o2i(o1(x)))

def test_suite_all():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_ops("test_linop"))
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite_all()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
