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

from ..import base
from ...patch import jsonpickle
from ...materials import compoundfromname
from ...materials import element
from ...materials import mixture


class test_base(unittest.TestCase):

    @unittest.skipIf(compoundfromname.xraylib is None,
                     "xraylib not installed")
    def test_serialize(self):
        ca = element.Element("Ca")
        calcite = compoundfromname.compoundfromname("calcite")
        hematite = compoundfromname.compoundfromname("hematite")
        goethite = compoundfromname.compoundfromname("goethite")
        mix = mixture.Mixture([ca, hematite, goethite, calcite],
                              [0.25, 0.25, 0.25, 0.25],
                              name="iron oxides")
        attenuators = {'A': ca, 'B': calcite, 'C': mix}

        d1 = base.Material()
        d2 = jsonpickle.loads(jsonpickle.dumps(d1))
        self.assertEqual(d1, d2)
        d1 = base.Material(attenuators=attenuators)
        d2 = jsonpickle.loads(jsonpickle.dumps(d1))
        self.assertEqual(d1, d2)

        d1 = base.SolidState()
        d2 = jsonpickle.loads(jsonpickle.dumps(d1))
        self.assertEqual(d1, d2)
        d1 = base.SolidState(attenuators=attenuators,
                             ehole=3.6)
        d2 = jsonpickle.loads(jsonpickle.dumps(d1))
        self.assertEqual(d1, d2)

        d1 = base.CentricCone()
        d2 = jsonpickle.loads(jsonpickle.dumps(d1))
        self.assertEqual(d1, d2)
        d1 = base.CentricCone(attenuators=attenuators,
                              ehole=3.6,
                              activearea=0.8)
        d2 = jsonpickle.loads(jsonpickle.dumps(d1))
        self.assertEqual(d1, d2)


def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_base("test_serialize"))
    return testSuite


if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
