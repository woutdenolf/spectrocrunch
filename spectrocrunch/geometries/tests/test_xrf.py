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

from ..import xrf
from ...utils import units
from ...patch import jsonpickle


class test_xrf(unittest.TestCase):

    def test_distance(self):
        geometry = xrf.factory("LinearXRFGeometry",
                               zerodistance=60.,
                               detectorposition=-10,
                               positionunits="mm",
                               detector=None,
                               source=None)

        self.assertEqual(geometry.distance.to("cm").magnitude, 5)
        self.assertEqual(geometry.detectorposition.to("cm").magnitude, -1)

        geometry.detectorposition = 30
        self.assertEqual(geometry.distance.to("cm").magnitude, 9)
        self.assertEqual(geometry.detectorposition.to("cm").magnitude, 3)

        geometry.detectorposition = units.Quantity(4, "cm")
        self.assertEqual(geometry.distance.to("cm").magnitude, 10)
        self.assertEqual(geometry.detectorposition.to("cm").magnitude, 4)

        geometry.distance = 4
        self.assertEqual(geometry.distance.to("cm").magnitude, 4)
        self.assertEqual(geometry.detectorposition.to("cm").magnitude, -2)

        geometry.calibrate_manually(10)
        self.assertEqual(geometry.distance.to("cm").magnitude, 10)
        self.assertEqual(geometry.detectorposition.to("cm").magnitude, -2)

        geometry.calibrate_manually(units.Quantity(20, "mm"))
        self.assertEqual(geometry.distance.to("cm").magnitude, 2)
        self.assertEqual(geometry.detectorposition.to("cm").magnitude, -2)

        geometry.distance = units.Quantity(50, "mm")
        self.assertEqual(geometry.distance.to("cm").magnitude, 5)
        self.assertEqual(geometry.detectorposition.to("cm").magnitude, 1)

    @unittest.skipIf(xrf.compoundfromname.xraylib is None,
                     "xraylib not installed")
    def test_serialize(self):
        exclude = 'XRFGeometry', 'LinearXRFGeometry'
        for name, cls in xrf.XRFGeometry.clsregistry.items():
            if name not in exclude:
                g1 = cls()
                g2 = jsonpickle.loads(jsonpickle.dumps(g1))
                self.assertEqual(g1, g2)


def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_xrf("test_distance"))
    testSuite.addTest(test_xrf("test_serialize"))
    return testSuite


if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
