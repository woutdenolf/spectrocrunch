# -*- coding: utf-8 -*-

import unittest

from ...detectors import area
from ...materials import multilayer
from ...geometries import flatarea
from ...sources import xray as xraysources
from ...detectors import diode
from ...materials import scintillators
from ...materials import lenses
from ...sources import emspectrum
from ...math import noisepropagation
from ...materials import compoundfromname
from ...utils.instance import isarray
from ...patch.pint import ureg

import numpy as np


class test_objects(unittest.TestCase):
    def _assertRV(self, RV):
        if isarray(RV):
            tmp = noisepropagation.E(RV)
            self.assertTrue(np.all(tmp == tmp[0]))
            tmp = noisepropagation.S(RV)
            self.assertTrue(np.all(tmp == tmp[0]))

    def _checkprop(self, o, **kwargs):
        self._checkprop1(o, **kwargs)
        self._checkprop2(o, **kwargs)

    def _checkprop1(self, o, vis=False, **kwargs):
        # Both arrays
        if vis:
            s = emspectrum.Discrete(ureg.Quantity([500, 500], "nm"), [1, 1])
        else:
            s = np.asarray([7.0, 7.0])
        N = noisepropagation.poisson([1e5, 1e5, 1e5])

        Nout = o.propagate(N, s, **kwargs)

        if vis:
            self.assertEqual(Nout.shape, (1, 3))
        else:
            self.assertEqual(Nout.shape, (2, 3))
        self._assertRV(Nout)

        # N array
        if vis:
            s = emspectrum.Discrete(ureg.Quantity(500, "nm"), 1)
        else:
            s = 7.0
        N = noisepropagation.poisson([1e5, 1e5, 1e5])

        Nout = o.propagate(N, s, **kwargs)

        self.assertEqual(Nout.shape, (1, 3))
        self._assertRV(Nout)

        # Energy array
        if vis:
            s = emspectrum.Discrete(ureg.Quantity([500, 500], "nm"), [1, 1])
        else:
            s = np.asarray([7.0, 7.0])
        N = noisepropagation.poisson(1e5)

        Nout = o.propagate(N, s, **kwargs)

        if vis:
            self.assertFalse(isarray(Nout))
        else:
            self.assertEqual(Nout.shape, (2, 1))
        self._assertRV(Nout)

        # Not arrays
        if vis:
            s = emspectrum.Discrete(ureg.Quantity(500, "nm"), 1)
        else:
            s = 7.0
        N = noisepropagation.poisson(1e5)

        Nout = o.propagate(N, s, **kwargs)

        self.assertFalse(isarray(Nout))
        self._assertRV(Nout)

    def _checkprop2(self, o, vis=False, **kwargs):
        if vis:
            s = emspectrum.Discrete(ureg.Quantity(500, "nm"), 1)
        else:
            s = 7.0
        N = noisepropagation.poisson(1e5)
        Nout = o.propagate(N, s, forward=True, **kwargs)
        N2 = o.propagate(Nout, s, forward=False, **kwargs)
        Nout2 = o.propagate(noisepropagation.E(N), s, forward=True, **kwargs)
        N3 = o.propagate(Nout2, s, forward=False, **kwargs)
        self.assertAlmostEqual(noisepropagation.E(N), noisepropagation.E(N2))
        self.assertAlmostEqual(noisepropagation.S(N), noisepropagation.S(N2))
        self.assertAlmostEqual(noisepropagation.E(N), N3)
        self.assertAlmostEqual(noisepropagation.E(Nout), Nout2)

    def test_detectors(self):
        self.assertRaises(RuntimeError, area.factory, "noclassname")

        src = xraysources.factory("synchrotron")
        detector = area.factory("PCO Edge 5.5")
        geometry = flatarea.factory("perpendicular", detector=detector, source=src)

        self._checkprop(detector, tframe=2, nframe=10, vis=True)

    @unittest.skipIf(lenses.visirlib.PyTMM is None, "PyTMM is not installed")
    def test_lenses(self):
        self.assertRaises(RuntimeError, lenses.factory, "noclassname")
        o = lenses.factory("Mitutoyo ID21 10x")
        self._checkprop(o, nrefrac=1.1, vis=True)
        o = lenses.factory("Mitutoyo ID21 10x")
        self._checkprop(o, nrefrac=1.1, vis=True)

    @unittest.skipIf(compoundfromname.xraylib is None, "xraylib is not installed")
    def test_scintillators(self):
        self.assertRaises(RuntimeError, scintillators.factory, "noclassname")
        o = scintillators.factory("GGG ID21", thickness=13)
        self._checkprop(o)

        # SNR = noisepropagation.SNR(N)

        # import matplotlib.pyplot as plt
        # plt.figure()
        # energy = np.linspace(2,9,100,dtype=float)
        # plt.plot(energy,o.transmission(energy))

        o = scintillators.factory("LSO ID21", thickness=10)
        self._checkprop(o)
        self._checkprop2(o)

        # plt.plot(energy,o.transmission(energy))
        # plt.show()

    @unittest.skipIf(compoundfromname.xraylib is None, "xraylib not installed")
    def test_materials(self):
        self.assertRaises(RuntimeError, multilayer.factory, "noclassname")

        detector = area.factory("PCO Edge 5.5")
        src = xraysources.factory("synchrotron")
        geometry = flatarea.factory("perpendicular", detector=detector, source=src)

        o = multilayer.Multilayer(
            material=compoundfromname.compoundfromname("ultralene"),
            thickness=4e-4,
            geometry=geometry,
        )
        self._checkprop(o)

    def test_diodes(self):
        return  # TODO
        o = diode.factory("idet")
        self._checkprop(o)


def main_test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_objects("test_materials"))
    testSuite.addTest(test_objects("test_scintillators"))
    testSuite.addTest(test_objects("test_detectors"))
    testSuite.addTest(test_objects("test_lenses"))
    testSuite.addTest(test_objects("test_diodes"))

    return testSuite


if __name__ == "__main__":
    import sys

    mysuite = main_test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
