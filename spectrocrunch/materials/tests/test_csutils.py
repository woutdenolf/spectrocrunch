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

import unittest

from ..import csutils
import numpy as np


class test_csutils(unittest.TestCase):

    @unittest.skipIf(csutils.xraylib.XRayInit is None,
                     "xraylib not installed")
    def test_shape(self):
        keep = csutils.xraylib.xraylib_np
        csutils.xraylib.xraylib_np = None
        try:
            self.assert_shape()
        finally:
            csutils.xraylib.xraylib_np = keep

    @unittest.skipIf(csutils.xraylib.xraylib_np is None,
                     "xraylib_np not installed")
    def test_shape_np(self):
        self.assert_shape()

    def assert_shape(self):
        method = 'CS_Total_Kissel'
        Z, E = 1, 7
        result = csutils.eval(method, Z, E)
        self.assertTrue(isinstance(result, float))
        Z, E = 1, np.linspace(7, 7.5, 5)
        result = csutils.eval(method, Z, E)
        self.assertEqual(result.shape, (5, ))
        result, postfunc = csutils.eval(method, Z, E, applypost=False)
        self.assertEqual(result.shape, (1, 5))


def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_csutils("test_shape"))
    testSuite.addTest(test_csutils("test_shape_np"))
    return testSuite


if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
