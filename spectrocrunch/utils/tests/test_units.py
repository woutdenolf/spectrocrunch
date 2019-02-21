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
import itertools
import numpy as np

from .gendata import gendata
from .. import units
from .. import instance
from ...patch.pint import ureg


class test_units(unittest.TestCase):

    def test_quantunits(self):
        for q, arr, expand in itertools.product([True, False], repeat=3):
            # x: scalar or array
            if arr:
                x = np.arange(1., 10)
            else:
                x = 10.

            # y: scalar, array, quantity(scalar), array(quantity), quantity(array)
            if q:
                unit0 = ureg.millimeter
                if arr and expand:
                    y = np.vectorize(lambda a: units.Quantity(
                        a, units=unit0), otypes=[object])(x)
                else:
                    y = units.Quantity(x, units=unit0)
                unit1 = ureg.meter
                munit = 1e-3
            else:
                unit0 = None
                unit1 = None
                munit = 1
                y = x

            # Test bubbling up and down
            if q:
                # z: quantity(scalar), quantity(array)
                z = units.Quantity(x, units=unit0)
                np.testing.assert_array_equal(
                    units.asqarray(y), units.asqarray(z))
                # z: array(quantity)
                z = np.vectorize(lambda a: units.Quantity(
                    a, units=unit0), otypes=[object])(x)
                np.testing.assert_array_equal(
                    units.Quantity(y), units.Quantity(z))

            # Test magnitude
            a = x*munit
            b = units.magnitude(y, unit1)
            if arr:
                np.testing.assert_array_equal(a, b)
            else:
                self.assertEqual(a, b)

            # Test unit conversion
            if q:
                a = units.Quantity(x, units=unit0)
                b = units.quantity_like(y, a)
            if arr:
                np.testing.assert_array_equal(a, b)
            else:
                self.assertEqual(a, b)

    def test_asqarray(self):
        from array import array
        for k, v in gendata().items():
            if instance.isarray(v) or instance.isscalar(v):
                varr = units.asqarray(v)
                msg = 'asqarray({}) is not an array'.format(k)
                self.assertTrue(instance.isqarray(varr), msg=msg)
        for k, v in gendata().items():
            if instance.isarray(v) or instance.isscalar(v):
                barr = instance.isarray(v)
                varr, restore = units.asqarrayf(v)
                msg = 'asqarrayf({}) is not an array'.format(k)
                self.assertTrue(instance.isqarray(varr), msg=msg)
                v2 = restore(varr)
                barr2 = instance.isqarray(v2)
                if barr:
                    msg = 'restore({}) is not an array'.format(k)
                else:
                    msg = 'restore({}) is an array'.format(k)
                self.assertEqual(barr, barr2, msg=msg)

def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_units("test_quantunits"))
    testSuite.addTest(test_units("test_asqarray"))
    return testSuite


if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
