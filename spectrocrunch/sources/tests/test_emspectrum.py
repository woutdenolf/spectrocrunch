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
import numpy as np

from .. import emspectrum
from ...patch.pint import ureg
from ...patch import jsonpickle


class test_emspectrum(unittest.TestCase):

    def test_discrete(self):
        s1 = emspectrum.Discrete(ureg.Quantity(
            [100, 200, 300], 'nm'), [1, 1, 1])
        s2 = emspectrum.Discrete(ureg.Quantity([200, 400], 'nm'), [1, 1])
        s3 = s1+s2
        def func(x): return ureg.Quantity(x, 'nm').to('keV', 'spectroscopy')
        np.testing.assert_array_equal(s3.energies, func([100, 200, 300, 400]))
        np.testing.assert_array_equal(s3.ratios, [0.2, 0.4, 0.2, 0.2])
        s4 = emspectrum.Discrete(ureg.Quantity([150, 350], 'nm'), [1, 2])
        self.assertEqual(s4.sample(s3), 1.5*1/3. + 1*2/3.)

    def test_dark(self):
        s1 = emspectrum.Dark()
        s4 = emspectrum.Discrete(ureg.Quantity([150, 350], 'nm'), [1, 2])
        self.assertEqual(s4.sample(s1), 0)

    def test_serialize(self):
        s1 = emspectrum.Discrete(ureg.Quantity([100, 300], 'nm'), [1, 1])
        s2 = jsonpickle.loads(jsonpickle.dumps(s1))
        self.assertEqual(s1, s2)
        s1 = emspectrum.Dark()
        s2 = jsonpickle.loads(jsonpickle.dumps(s1))
        self.assertEqual(s1, s2)


def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_emspectrum("test_discrete"))
    testSuite.addTest(test_emspectrum("test_dark"))
    testSuite.addTest(test_emspectrum("test_serialize"))
    return testSuite


if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
