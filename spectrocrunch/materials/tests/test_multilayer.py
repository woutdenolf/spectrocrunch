# -*- coding: utf-8 -*-

import unittest
import numpy as np
import scipy.integrate
import scipy.special

from .. import multilayer
from .. import compoundfromformula
from .. import compoundfromname
from .. import mixture
from .. import element
from .. import types
from .. import xrayspectrum
from ...geometries import xrf as xrfgeometries
from ...detectors import xrf as xrfdetectors
from ...sources import xray as xraysources
from ...utils import timing
from ...utils import listtools
from ...patch import jsonpickle


class test_multilayer(unittest.TestCase):

    def test_expi(self):
        for x in [1e-4, 1, 10, np.inf]:  # x > 0
            int1 = scipy.integrate.quad(
                lambda a: np.exp(-x/np.cos(a))*np.tan(a), 0, np.pi/2)
            int2 = -scipy.special.expi(-x)
            int3 = scipy.special.exp1(x)
            self.assertGreaterEqual(int2, int1[0]-int1[1])
            self.assertLessEqual(int2, int1[0]+int1[1])
            self.assertGreaterEqual(int3, int1[0]-int1[1])
            self.assertLessEqual(int3, int1[0]+int1[1])

    def _multilayer1(self):
        src = xraysources.factory("synchrotron")
        detector = xrfdetectors.factory("XRFDetector", activearea=0.50)
        geometry = xrfgeometries.factory("sxm120",
                                         detectorposition=-15.,
                                         detector=detector,
                                         source=src)
        geometry.solidangle = 4*np.pi*0.1
        c1 = compoundfromformula.CompoundFromFormula(
            "CaCO3", 2.71, name="calcite")
        c2 = compoundfromformula.CompoundFromFormula(
            "Fe2O3", 5.3, name="hematite")
        c3 = compoundfromformula.CompoundFromFormula(
            "PbCO3", 6.53, name="cerussite")
        l = [c1, mixture.Mixture([c2, c3], [1, 1], types.fraction.mole)]
        thickness = [10e-4, 20e-4]
        o = multilayer.Multilayer(
            material=l, thickness=thickness, geometry=geometry)
        return o

    def _multilayer2(self):
        element1 = element.Element("Ca")
        compound1 = compoundfromformula.CompoundFromFormula(
            "PbSO4", density=6.29)
        compound2 = compoundfromformula.CompoundFromFormula(
            "CaSO4", density=2.32)
        mixture1 = mixture.Mixture([compound1, compound2], [
                                   0.5, 0.5], types.fraction.mass)
        src = xraysources.factory("synchrotron")
        detector = xrfdetectors.factory("leia")
        geometry = xrfgeometries.factory("sxm120",
                                         detectorposition=-15.,
                                         detector=detector,
                                         source=src)
        o = multilayer.Multilayer(material=[element1, compound1, mixture1],
                                  thickness=[1e-4, 1e-4, 1e-4],
                                  geometry=geometry)
        return o

    def _multilayer3(self):
        element1 = element.Element("Ca")
        src = xraysources.factory("synchrotron")
        detector = xrfdetectors.factory("leia")
        geometry = xrfgeometries.factory("sxm120",
                                         detectorposition=4.,
                                         detector=detector,
                                         source=src)
        o = multilayer.Multilayer(material=[element1],
                                  thickness=[1e-4],
                                  geometry=geometry)
        return o

    @unittest.skipIf(xrfdetectors.compoundfromname.xraylib is None, "xraylib not installed")
    def test_attenuationinfo(self):
        src = xraysources.factory("synchrotron")
        detector = xrfdetectors.factory("leia")
        geometry = xrfgeometries.factory("linearxrfgeometry", anglein=90, angleout=45,
                                         azimuth=0, detectorposition=0., positionunits="mm", detector=detector, source=src)
        o = multilayer.Multilayer(material=[compoundfromname.compoundfromname("hematite"),
                                            compoundfromname.compoundfromname(
                                                "hydrocerussite"),
                                            compoundfromname.compoundfromname("calcite")],
                                  thickness=[9e-4, 9e-4, 9e-4], geometry=geometry)

        with o.cachectx("layerinfo"):

            for energy in [8, np.array([8]), np.array([8, 8.5])]:

                with o.cachectx("attenuationinfo", energy):
                    nenergy = np.asarray(energy).size

                    rho = o.density[:, np.newaxis]
                    d = o.xraythicknessin[:, np.newaxis]
                    mu = o.mass_att_coeff(energy)
                    murhod = mu*rho*d

                    nsub = 3
                    n = nsub*o.nlayers+1
                    A = np.zeros((n, nenergy))
                    np.cumsum(np.repeat(murhod, nsub, axis=0) /
                              nsub, out=A[1:, :], axis=0)

                    z = np.zeros(n)
                    np.cumsum(np.repeat(d, nsub)/nsub, out=z[1:])

                    np.testing.assert_allclose(np.squeeze(
                        A), np.squeeze(o._cum_attenuation(z, energy)))

                    for i in range(n):
                        np.testing.assert_allclose(A[i, :], np.squeeze(
                            o._cum_attenuation(z[i], energy)))
                        cosaij = np.ones(n)
                        cosaij[0:i] = -1

                        T1 = np.exp(-np.abs(A-A[i, np.newaxis, :]))
                        T2 = o._transmission(z[i], z, cosaij, energy)
                        np.testing.assert_allclose(
                            np.squeeze(T1), np.squeeze(T2))

                        for j in range(n):
                            T1 = np.exp(-abs(A[j, :]-A[i, :]))
                            T2 = np.squeeze(o._transmission(
                                z[i], z[j], cosaij[j], energy))
                            np.testing.assert_allclose(T1, T2)

    def assertSpectrumEqual(self, spectrum1, spectrum2, rtol=1e-07, compfisx=False):
        if compfisx:
            # Fisx does not return scattering
            for k in ["Rayleigh", "Compton"]:
                try:
                    spectrum1.pop(k)
                except KeyError:
                    pass
                try:
                    spectrum2.pop(k)
                except KeyError:
                    pass
            # Fisx adds sources
            for k in spectrum1:
                spectrum1[k] = np.sum(spectrum1[k])
                spectrum2[k] = np.sum(spectrum2[k])
        self.assertEqual(set(spectrum1.keys()), set(spectrum2.keys()))
        for k in spectrum1:
            np.testing.assert_allclose(spectrum1[k], spectrum2[k], rtol=rtol)

    @unittest.skipIf(xrfdetectors.compoundfromname.xraylib is None, "xraylib not installed")
    def test_primary_simple(self):
        # Define setup: single layer
        compound1 = compoundfromformula.CompoundFromFormula(
            "Mn1Fe2Ca3O4", density=7.)
        src = xraysources.factory("synchrotron")
        detector = xrfdetectors.factory("leia")
        geometry = xrfgeometries.factory(
            "sxm120", detectorposition=-15., detector=detector, source=src)
        detector.attenuators["BeamFilter"] = {
            "material": element.Element('Si'), "thickness": 10e-4}
        detector.attenuators["Filter"] = {
            "material": element.Element('Si'), "thickness": 10e-4}
        o = multilayer.Multilayer(material=[compound1], thickness=[
                                  50e-4], geometry=geometry)

        integratormult = geometry.solidangle/geometry.cosnormin  # srad
        for anglesign in [1, -1]:
            detector.geometry.angleout = anglesign*abs(geometry.angleout)
            energy0 = [8, 9]
            weights = [1, 2]

            # Calculate spectrum manually (1 layer, 1 interaction)
            emax = max(energy0)
            rates = dict(compound1.xrayspectrum(
                energy0, emin=2, emax=emax).probabilities)  # 1/cm
            mu0 = compound1.mass_att_coeff(energy0)
            for k in rates:
                energy1 = k.energy(
                    **detector.geometry.xrayspectrumkwargs())
                mu1 = compound1.mass_att_coeff(energy1)
                chi = mu1/o.geometry.cosnormout-mu0/o.geometry.cosnormin
                chi = chi*compound1.density
                thicknesscor = (np.exp(chi*o.thickness[0])-1)/chi  # cm
                if not geometry.reflection:
                    thicknesscor *= np.exp(-mu1*compound1.density *
                                           o.thickness[0]/o.geometry.cosnormout)
                eff = o.geometry.efficiency(energy0, energy1)
                if not isinstance(k, xrayspectrum.FluoZLine):
                    eff = np.diag(eff)
                eff = np.squeeze(eff)
                rates[k] = rates[k]*eff*integratormult * \
                    thicknesscor*weights/np.sum(weights)
            spectrumRef = xrayspectrum.Spectrum()
            spectrumRef.update(rates)
            spectrumRef.xlim = [2, emax]
            spectrumRef.density = 1
            spectrumRef.title = str(self)
            spectrumRef.type = spectrumRef.TYPES.rate

            # Compare with different methods
            for asarray in False, True:
                if asarray:
                    energy0 = [energy0]*2
                    weights = [weights]*2
                with timing.timeit_logger(name="analytical"):
                    spectra1 = o.xrayspectrum(
                        energy0, weights=weights, emin=2, emax=emax,
                        method="analytical", ninteractions=1)
                with timing.timeit_logger(name="numerical"):
                    spectra2 = o.xrayspectrum(
                        energy0, weights=weights, emin=2, emax=emax,
                        method="numerical", ninteractions=1)
                with timing.timeit_logger(name="fisx"):
                    spectra3 = o.xrayspectrum(
                        energy0, weights=weights, emin=2, emax=emax,
                        method="fisx", ninteractions=1)
                if not asarray:
                    spectra1 = [spectra1]
                    spectra2 = [spectra2]
                    spectra3 = [spectra3]
                for spectrum1 in spectra1:
                    self.assertSpectrumEqual(
                        spectra1[0], spectrum1, rtol=1e-07)
                for spectrum1 in spectra1:
                    self.assertSpectrumEqual(
                        spectrum1, spectrumRef, rtol=1e-06)
                    for spectrum2 in spectra2:
                        self.assertSpectrumEqual(
                            spectrum1, spectrum2, rtol=1e-07)
                    for spectrum3 in spectra3:
                        self.assertSpectrumEqual(
                            spectrum1, spectrum3, rtol=1e-06, compfisx=True)

    @unittest.skipIf(xrfdetectors.compoundfromname.xraylib is None,
                     "xraylib not installed")
    def test_primary_complex(self):
        compound1 = compoundfromformula.CompoundFromFormula("MnFe", density=6.)
        compound2 = compoundfromformula.CompoundFromFormula(
            "Mn1Fe2Ca3O4", density=9.)
        src = xraysources.factory("synchrotron")
        detector = xrfdetectors.factory("leia")
        geometry = xrfgeometries.factory(
            "sxm120", detectorposition=-15., detector=detector, source=src)
        detector.attenuators["BeamFilter"] = {
            "material": element.Element('Si'), "thickness": 10e-4}
        detector.attenuators["Filter"] = {
            "material": element.Element('Si'), "thickness": 10e-4}
        o = multilayer.Multilayer(material=[compound1, compound2],
                                  thickness=[5e-4, 10e-4],
                                  geometry=geometry)

        for anglesign in [1, -1]:
            detector.geometry.angleout = anglesign*abs(geometry.angleout)
            energy0 = [8, 9]
            weights = [1, 2]
            with timing.timeit_logger(name="analytical"):
                spectrum1 = o.xrayspectrum(
                    energy0, weights=weights, emin=2,
                    method="analytical", ninteractions=1)
            with timing.timeit_logger(name="numerical"):
                spectrum2 = o.xrayspectrum(
                    energy0, weights=weights, emin=2,
                    method="numerical", ninteractions=1)
            with timing.timeit_logger(name="fisx"):
                spectrum3 = o.xrayspectrum(
                    energy0, weights=weights, emin=2,
                    method="fisx", ninteractions=1)
            self.assertSpectrumEqual(
                spectrum1, spectrum2, rtol=1e-04)
            self.assertSpectrumEqual(
                spectrum1, spectrum3, rtol=1e-04, compfisx=True)

    @unittest.skip("TODO")
    def test_secondary(self):
        compound1 = compoundfromformula.CompoundFromFormula("MnFe", density=6.)

        compound2 = compoundfromformula.CompoundFromFormula(
            "Mn1Fe2Ca3O4", density=9.)

        src = xraysources.factory("synchrotron")
        detector = xrfdetectors.factory("leia")
        geometry = xrfgeometries.factory(
            "sxm120", detectorposition=-15., detector=detector, source=src)

        o = multilayer.Multilayer(material=[compound1, compound2],
                                  thickness=[1e-4, 1e-4],
                                  geometry=geometry)

        for anglesign in [-1]:
            detector.geometry.angleout = anglesign*abs(geometry.angleout)
            energy0 = [8, 9]
            weights = [1, 2]
            spectrum2 = o.xrayspectrum(
                energy0, weights=weights, emin=2, emax=energy0, method="numerical", ninteractions=2)
            spectrum3 = o.xrayspectrum(
                energy0, weights=weights, emin=2, emax=energy0, method="fisx", ninteractions=2)
            self.assertSpectrumEqual(
                spectrum2, spectrum3, rtol=1e-02, compfisx=True)

    @unittest.skipIf(xrfdetectors.compoundfromname.xraylib is None,
                     "xraylib not installed")
    def test_serialize(self):
        m1 = self._multilayer2()
        m2 = jsonpickle.loads(jsonpickle.dumps(m1))
        self.assertEqual(m1, m2)

    @unittest.skipIf(xrfdetectors.compoundfromname.xraylib is None,
                     "xraylib not installed")
    def test_mixlayers(self):
        m = self._multilayer1()
        energy = 10
        t = sum(m.xraythicknessin)
        mix = m.mixlayers()
        abs1 = m.absorbance(energy)
        abs2 = mix.mass_att_coeff(energy)*mix.density*t
        self.assertAlmostEqual(abs1, abs2)


def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_multilayer("test_expi"))
    testSuite.addTest(test_multilayer("test_attenuationinfo"))
    testSuite.addTest(test_multilayer("test_primary_simple"))
    testSuite.addTest(test_multilayer("test_primary_complex"))
    testSuite.addTest(test_multilayer("test_serialize"))
    testSuite.addTest(test_multilayer("test_secondary"))
    testSuite.addTest(test_multilayer("test_mixlayers"))
    return testSuite


if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
