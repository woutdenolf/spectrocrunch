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
import numpy as np

from ..import lut
from .. import units
from ...patch import jsonpickle


class test_lut(unittest.TestCase):

    def test_sort(self):
        l1 = lut.LUT()
        l1.add([3,2],[0,-1])
        l1.add([1,0],[-2,-3])
        self.assertEqual(l1.x.magnitude.tolist(), list(range(4)))
        self.assertEqual(l1.y.magnitude.tolist(), list(range(-3,1)))

    def test_units(self):
        l1 = lut.LUT()
        l1.add(units.Quantity([3,2], 'cm'),
               units.Quantity([0,-1], 'keV'))
        l1.add(units.Quantity([10,0], 'mm'),
               units.Quantity([-2000,-3000], 'eV'))
        self.assertEqual(l1.x.magnitude.tolist(), list(range(4)))
        self.assertEqual(l1.y.magnitude.tolist(), list(range(-3,1)))
        self.assertEqual(l1.x.units, units.ureg.Unit('cm'))
        self.assertEqual(l1.y.units, units.ureg.Unit('keV'))

    def test_interpolate(self):
        l1 = lut.LUT(kind='linear')
        l1.add(units.Quantity(7,'keV'), units.Quantity(10,'mm'))
        l1.add(units.Quantity(7400,'eV'), units.Quantity(2,'cm'))
        func = lambda x: l1(units.Quantity(x,'keV')).to('mm').magnitude
        self.assertEqual(func(7.2), 15)
        np.testing.assert_allclose(func([7.2]), [15])
        np.testing.assert_allclose(func([7.1, 7.2, 7.3]),
                                   [12.5, 15, 17.5])

    def test_serialize(self):
        l1 = lut.LUT()
        l2 = jsonpickle.decode(jsonpickle.encode(l1))
        self.assertEqual(l1, l2)
        l1.add(1, 2)
        l2 = jsonpickle.decode(jsonpickle.encode(l1))
        self.assertEqual(l1, l2)
        l1.add([4, 5], [6, 7])
        l2 = jsonpickle.decode(jsonpickle.encode(l1))
        self.assertEqual(l1, l2)


def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_lut("test_sort"))
    testSuite.addTest(test_lut("test_units"))
    testSuite.addTest(test_lut("test_interpolate"))
    testSuite.addTest(test_lut("test_serialize"))
    return testSuite


if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
