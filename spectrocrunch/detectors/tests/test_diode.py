# -*- coding: utf-8 -*-

import unittest
import numpy as np
from scipy import interpolate
from scipy import constants
import random
import itertools

from .. import diode
from ...resources import resource_filename
from ...optics import xray as xrayoptics
from ...patch.pint import ureg
from ...patch import jsonpickle


class test_diode(unittest.TestCase):

    def _spec_idet_cpstoflux(self, ph_E, ph_I, ph_gain):
        """
        Photon calculation in spec

        Args:
            ph_E(array): keV
            ph_I(array): idet counts
            ph_gain(array): V/A
        """
        PH_DIST = 0
        ph_coeff = 0
        ph_I = ph_I * 10.**(-5-ph_gain)
        ptb = np.loadtxt(resource_filename('id21/ptb.dat'))
        fptb = interpolate.interp1d(ptb[:, 0], ptb[:, 1])
        ird = np.loadtxt(resource_filename('id21/ird.dat'))
        fird = interpolate.interp1d(ird[:, 0], ird[:, 1])
        ph_PTB = fptb(ph_E)
        ph_factor = fird(ph_E)
        ph_calib = ph_factor * ph_PTB
        return ph_I / (ph_E * constants.elementary_charge * np.exp(-ph_coeff * PH_DIST) * ph_calib)

    @unittest.skipIf(diode.compoundfromname.xraylib is None, "xraylib not installed")
    def test_idet(self):
        o1 = diode.SXM_IDET(model=True)
        o2 = diode.SXM_IDET(model=False)
        self.assertAlmostEqual(
            o1._chargepersamplephoton(5.2).to('coulomb').magnitude,
            o2._chargepersamplephoton(5.2).to('coulomb').magnitude)

        for model in [True, False]:
            o = diode.SXM_IDET(model=model)
            o.gain = 1e8  #ureg.Quantity(1e8, 'V/A')
            self._assertDiode(o)

            # Compare to LUT SXM calculation
            energy = np.arange(3, 7)[:, np.newaxis]  # several single-energy sources
            I = np.arange(5, 10)[np.newaxis, :]*ureg.Quantity(1e5, "Hz")  # several fluxes for each sour
            np.testing.assert_array_almost_equal(
                o.fluxtocps(energy, o.cpstoflux(energy, I)),
                np.repeat(I, energy.size, axis=0))
            flux1 = self._spec_idet_cpstoflux(energy, I.magnitude, 8)
            flux2 = o.cpstoflux(energy, I)
            if model:
                for f1, f2 in zip(flux1.flatten(), flux2.flatten()):
                    np.testing.assert_approx_equal(
                        f1, f2.magnitude, significant=1)
            else:
                # 3% difference with spec
                np.testing.assert_allclose(flux1, flux2, rtol=0.03)

    def _assertDiode(self, o):
        energy = np.arange(3, 5)[:, np.newaxis]  # several single-energy sources
        cps = np.linspace(1e6, 1e7, 5)[np.newaxis, :]
        flux1 = o.cpstoflux(energy, cps)
        energy = np.arange(3, 6)[:, np.newaxis]  # several single-energy sources
        flux2 = o.cpstoflux(energy, cps)
        np.testing.assert_allclose(flux1, flux2[:-1, :])

        # Test inversions:
        o2 = o.op_cpstocurrent()*o.op_currenttocps()
        self.assertAlmostEqual(o2.m.magnitude, 1.)
        self.assertAlmostEqual(o2.b.magnitude, 0.)
        self.assertEqual(o2.m.units, ureg.dimensionless)
        self.assertEqual(o2.b.units, ureg.ampere)

        energy = np.arange(3, 5)[:, np.newaxis]
        o2 = o.op_fluxtocurrent(energy)*o.op_currenttoflux(energy)
        np.testing.assert_array_almost_equal(o2.m, 1.)
        np.testing.assert_array_almost_equal(o2.b.magnitude, 0.)
        self.assertEqual(o2.m.units, ureg.dimensionless)
        self.assertEqual(o2.b.units, ureg.ampere)

        o2 = o.op_fluxtocps(energy)*o.op_cpstoflux(energy)
        np.testing.assert_array_almost_equal(o2.m, 1.)
        np.testing.assert_array_almost_equal(o2.b.magnitude, 0.)
        self.assertEqual(o2.m.units, ureg.dimensionless)
        self.assertEqual(o2.b.units, ureg.hertz)

    @unittest.skipIf(diode.compoundfromname.xraylib is None,
                     "xraylib not installed")
    def test_iodet1(self):
        parameters = [(True, False,),
                      ("solidangle", "thickness", "optics")]
        for params in itertools.product(*parameters):
            simplecalibration, caliboption = params
            if simplecalibration:
                # Only LUT changes when calibration
                if caliboption != "optics":
                    continue
            # First optic with LUT will be used for calibration if caliboption == "optics"
            airpath = xrayoptics.Filter(material="air", thickness=10)
            optics = [airpath, "kb"]
            o = diode.SXM_IODET1(optics=optics,
                                 simplecalibration=simplecalibration)
            o.thickness = 1e-4
            o.geometry.solidangle = 0.1
            if caliboption == "solidangle":
                # Same counts but the flux is * m -> diode yield (solid angle) is / m
                energy = 7.
                cps = np.linspace(1e4, 1e5, 10)
                sampleflux = o.cpstoflux(energy, cps)
                m = 2.
                sa = o.geometry.solidangle
                flux2 = sampleflux * m
                o.calibrate(cps, flux2, energy, caliboption=caliboption)
                flux1 = o.cpstoflux(energy, cps)
                np.testing.assert_allclose(flux1, flux2)
                np.testing.assert_allclose(o.geometry.solidangle, sa/m)
            elif caliboption == "optics":
                # Same counts but the flux is * m -> optics transmission is * m
                energy = [[7.], [7.1], [7.2]]
                cps = np.linspace(1e4, 1e5, 10)
                sampleflux = o.cpstoflux(energy, cps[np.newaxis, :])

                m = np.array([[0.5], [(0.5 + 0.2)/2.], [0.2]])
                flux2 = sampleflux * m
                o.calibrate(cps, flux2[0], energy[0],
                            caliboption=caliboption)
                o.calibrate(cps, flux2[-1], energy[-1],
                            caliboption=caliboption)
                flux1 = o.cpstoflux(energy, cps)
                if simplecalibration:
                    np.testing.assert_allclose(flux1[0], flux2[0])
                    # TODO: should this be the case?
                    #np.testing.assert_allclose(cps, o.fluxtocps(energy[1], flux2[1]))
                    np.testing.assert_allclose(flux1[-1], flux2[-1])
                else:
                    np.testing.assert_allclose(flux1, flux2)
                    np.testing.assert_allclose(
                        o.caliboptic.transmission(energy), m)
            elif caliboption == "thickness":
                energy = 7.
                cps = np.linspace(1e4, 1e5, 10)
                sampleflux = o.cpstoflux(energy, cps)
                m = 2.
                flux2 = sampleflux * m
                o.calibrate(cps, flux2, energy, caliboption=caliboption)
                np.testing.assert_allclose(flux1, flux2)
            self._assertDiode(o)

    @unittest.skipIf(diode.compoundfromname.xraylib is None,
                     "xraylib not installed")
    def test_serialize(self):
        exclude = 'PNdiode', 'CalibratedPNdiode', 'NonCalibratedPNdiode', 'SXM_PTB'
        for name, cls in diode.PNdiode.clsregistry.items():
            if name not in exclude:
                d1 = cls()
                d2 = jsonpickle.loads(jsonpickle.dumps(d1))
                self.assertEqual(d1, d2)


def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_diode("test_idet"))
    testSuite.addTest(test_diode("test_iodet1"))
    testSuite.addTest(test_diode("test_serialize"))
    return testSuite


if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
