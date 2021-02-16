# -*- coding: utf-8 -*-

import unittest

from .. import xrf
from ...utils import units
from ...patch import jsonpickle


class test_xrf(unittest.TestCase):
    def test_distance(self):
        geometry = xrf.factory(
            "LinearXRFGeometry",
            zerodistance=60.0,
            detectorposition=-10,
            positionunits="mm",
            detector=None,
            source=None,
        )

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

    @unittest.skipIf(xrf.compoundfromname.xraylib is None, "xraylib not installed")
    def test_serialize(self):
        exclude = "XRFGeometry", "LinearXRFGeometry"
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


if __name__ == "__main__":
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
