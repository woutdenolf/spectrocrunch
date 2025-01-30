# -*- coding: utf-8 -*-

import unittest
import numpy as np

from .. import qxrf
from ...utils import units
from ...patch import jsonpickle


class test_qxrf(unittest.TestCase):
    def geometryinstance(self):
        energy = 10
        geometryinstance = qxrf.factory("sxm1", simplecalibration=False)
        info = {
            "I0_counts": 300,
            "It_counts": 30,
            "time": 1,
            "dark": True,
            "gaindiodeI0": 1e8,
            "gaindiodeIt": 1e7,
        }
        geometryinstance.calibrate_diodes(**info)
        info["I0_counts"] = 10000
        info["It_counts"] = 100000
        info["energy"] = energy - 2
        info["dark"] = False
        geometryinstance.calibrate_diodes(**info)
        info["I0_counts"] = 5000
        info["energy"] = energy + 2
        geometryinstance.calibrate_diodes(**info)
        geometryinstance.reference = units.Quantity(1e9, "Hz")
        geometryinstance.defaultexpotime = units.Quantity(100, "ms")
        return geometryinstance

    @unittest.skipIf(
        qxrf.xrfdetectors.compoundfromname.xraylib is None, "xraylib not installed"
    )
    def test_flux(self):
        geometryinstance = self.geometryinstance()

        energy = 10
        time = 0.2
        refflux = 1e9

        flux = np.linspace(1e9, 1e8, 20)  # ph/sec
        iodet = geometryinstance.fluxtocps(energy, flux) * time

        flux2 = geometryinstance.responsetoflux(energy, iodet / time)
        np.testing.assert_allclose(flux, flux2)

        # Normalize data to the real flux (use flux reference)
        rates = np.random.poisson(np.full_like(flux, 100))  # 1/ph/sec
        data = flux * time * rates  # measured xrf
        # measured when flux whould have been refflux at each poi1e9
        dataref = refflux * time * rates
        ref = units.Quantity(refflux, "hertz")

        op, _, _, _ = geometryinstance.xrfnormop(energy, expotime=time, reference=ref)
        np.testing.assert_allclose(dataref, data / op(iodet))

        # Normalize data to the real flux (use iodet reference)
        iodetref = geometryinstance.fluxtocps(energy, refflux) * time
        ref = units.Quantity(iodetref, "dimensionless")

        op, _, _, _ = geometryinstance.xrfnormop(energy, expotime=time, reference=ref)
        np.testing.assert_allclose(dataref, data / op(iodet))

    @unittest.skipIf(
        qxrf.xrfgeometries.compoundfromname.xraylib is None, "xraylib not installed"
    )
    def test_serialize(self):
        xrfgeometries = []
        for detectorposition in [0, 1]:
            geometry = {
                "name": "sxm120",
                "parameters": {"detectorposition": detectorposition},
            }
            detector = {"name": "leia", "parameters": {}}
            xrfgeometries.append((geometry, detector))
        g1 = qxrf.factory("sxm1", xrfgeometries=xrfgeometries)
        g2 = jsonpickle.loads(jsonpickle.dumps(g1))
        self.assertEqual(g1, g2)

        exclude = ("QXRFGeometry",)
        for name, cls in qxrf.QXRFGeometry.clsregistry.items():
            if name not in exclude:
                g1 = cls()
                g2 = jsonpickle.loads(jsonpickle.dumps(g1))
                self.assertEqual(g1, g2)


def main_test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_qxrf("test_flux"))
    testSuite.addTest(test_qxrf("test_serialize"))
    return testSuite


if __name__ == "__main__":
    import sys

    mysuite = main_test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
