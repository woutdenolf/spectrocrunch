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

from ..import xray
from ...utils import units
from ...patch import jsonpickle


class test_xray(unittest.TestCase):

    def test_interpolate(self):
        o1 = xray.KB(kind='linear')
        o1.set_transmission(7, 0.2)
        o1.set_transmission(units.Quantity(7400, 'eV'), 0.8)
        def transmission(x): return o1.transmission(units.Quantity(x, 'keV'))
        self.assertEqual(transmission(7.2), 0.5)
        np.testing.assert_allclose(transmission([7.2]), [0.5])
        np.testing.assert_allclose(transmission([7.1, 7.2, 7.3]),
                                   [0.35, 0.5, 0.65])

    def test_serialize(self):
        exclude = ()
        for name, cls in xray.XrayOptics.clsregistry.items():
            if name not in exclude:
                o1 = cls()
                o2 = jsonpickle.decode(jsonpickle.encode(o1))
                self.assertEqual(o1, o2)


def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_xray("test_interpolate"))
    testSuite.addTest(test_xray("test_serialize"))
    return testSuite


if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
