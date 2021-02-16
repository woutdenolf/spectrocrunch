# -*- coding: utf-8 -*-

import unittest

from ...patch.pint import ureg
from ...patch import jsonpickle
from ...utils import instance
from .. import compoundfromformula
from .. import mixture
from .. import types
from ...geometries import xrf as xrfgeometries
from ...sources import xray as xraysources
from ...detectors import xrf as xrfdetectors

import numpy as np


class test_mixture(unittest.TestCase):
    @unittest.skipIf(
        xrfdetectors.compoundfromname.xraylib is None, "xraylib not installed"
    )
    def test_molefractions(self):
        compound1 = compoundfromformula.CompoundFromFormula("La2O3", 0)
        compound2 = compoundfromformula.CompoundFromFormula("SrO", 0)
        compound3 = compoundfromformula.CompoundFromFormula("Co2O3", 0)
        compound4 = compoundfromformula.CompoundFromFormula("Fe2O3", 0)

        mole = np.array([1.3, 2, 0.2, 3], dtype=float)
        m = mixture.Mixture(
            [compound1, compound2, compound3, compound4], mole, types.fraction.mole
        )

        # Test compound mole fractions
        nfrac1 = m.equivalents()
        nfrac2 = mole
        labels = ["La2O3", "SrO", "Co2O3", "Fe2O3"]
        for i in range(len(labels)):
            self.assertAlmostEqual(nfrac1[labels[i]], nfrac2[i])

        # Test elemental mole fractions
        nfrac1 = m.elemental_equivalents()
        nLa = 2 * mole[0]
        nSr = 1 * mole[1]
        nCo = 2 * mole[2]
        nFe = 2 * mole[3]
        nO = 3 * mole[0] + 1 * mole[1] + 3 * mole[2] + 3 * mole[3]
        nfrac2 = np.array([nSr, nCo, nFe, nO, nLa])
        labels = ["Sr", "Co", "Fe", "O", "La"]
        for i in range(len(labels)):
            self.assertAlmostEqual(nfrac1[labels[i]], nfrac2[i])

    @unittest.skipIf(
        xrfdetectors.compoundfromname.xraylib is None, "xraylib not installed"
    )
    def test_addcompounds(self):
        c1 = compoundfromformula.CompoundFromFormula("Co2O3", 1.5)
        c2 = compoundfromformula.CompoundFromFormula("Fe2O3", 1.6)
        m1 = mixture.Mixture([c1, c2], [2, 3], types.fraction.mole)
        n1 = m1.molefractions()
        snfrac1 = sum(m1.equivalents().values())
        for fractype in ["mass", "mole", "volume"]:
            m2 = mixture.Mixture([c1], [5], types.fraction.mole)
            if fractype == "mass":
                m2.addcompounds(c2, m1.massfractions()[c2], types.fraction.mass)
            elif fractype == "volume":
                m2.addcompounds(c2, m1.volumefractions()[c2], types.fraction.volume)
            else:
                m2.addcompounds(c2, n1[c2], types.fraction.mole)
            n2 = m2.molefractions()
            snfrac2 = sum(m2.equivalents().values())
            self.assertEqual(set(n1.keys()), set(n2.keys()))
            if fractype == "mole":
                self.assertEqual(snfrac1, snfrac2)
            else:
                self.assertAlmostEqual(1, snfrac2)
            self.assertAlmostEqual(m1.density, m2.density)
            for k in n1:
                self.assertAlmostEqual(n1[k], n2[k])

    @unittest.skipIf(
        xrfdetectors.compoundfromname.xraylib is None, "xraylib not installed"
    )
    def test_tocompound(self):
        c1 = compoundfromformula.CompoundFromFormula("Co2O3", 1.5)
        c2 = compoundfromformula.CompoundFromFormula("Fe2O3", 1.6)
        c3 = mixture.Mixture([c1, c2], [2, 3], types.fraction.mole).tocompound("mix")
        c4 = compoundfromformula.CompoundFromFormula("Co4O15Fe6", c3.density)

        n3 = c3.molefractions()
        n4 = c4.molefractions()
        self.assertEqual(sorted(n3.keys()), sorted(n4.keys()))
        for k in n3:
            self.assertEqual(n3[k], n4[k])

    def _spectrum(self, lstlines, thickness):
        src = xraysources.factory("synchrotron")
        detector = xrfdetectors.factory("leia")
        geometry = xrfgeometries.factory(
            "sxm120", detectorposition=-15.0, detector=detector, source=src
        )

        out = {}
        for spectrum, d in zip(lstlines, thickness):

            for line, prob in spectrum.probabilities:
                intensity = prob * detector.geometry.solidangle * d
                energy = line.energy(**detector.geometry.xrayspectrumkwargs())
                energy = instance.asarray(energy)
                intensity = instance.asarray(intensity)
                for en, inten in zip(energy, intensity):
                    if en in out:
                        out[en] += intensity
                    else:
                        out[en] = intensity
        return out

    def _spectrum_equal(self, spectrum1, spectrum2):
        np.testing.assert_allclose(sorted(spectrum1.keys()), sorted(spectrum2.keys()))

        for k in spectrum1:
            np.testing.assert_allclose(spectrum1[k], spectrum2[k])

    @unittest.skipIf(
        xrfdetectors.compoundfromname.xraylib is None, "xraylib not installed"
    )
    def test_cross_sections(self):
        c = [
            compoundfromformula.CompoundFromFormula("Co2O3", 1.5),
            compoundfromformula.CompoundFromFormula("Fe2O3", 1.6),
        ]
        rhoc = np.array([x.density for x in c])
        dc = np.array([0.1, 0.2])

        rhom = sum(rhoc * dc) / sum(dc)
        wc = rhoc * dc
        wc /= sum(wc)

        c3 = mixture.Mixture(c, wc, types.fraction.mass)
        c4 = c3.tocompound("mix")

        np.testing.assert_allclose(c3.density, rhom)
        np.testing.assert_allclose(c4.density, rhom)

        energy = 10.0
        muc = [x.mass_att_coeff(energy) for x in c]

        muL = sum(muc * rhoc * dc)
        np.testing.assert_allclose(
            c3.density * c3.mass_att_coeff(energy) * sum(dc), muL
        )
        np.testing.assert_allclose(
            c4.density * c4.mass_att_coeff(energy) * sum(dc), muL
        )

        spectrum = self._spectrum([x.xrayspectrum(energy) for x in c], dc)
        spectrum3 = self._spectrum([c3.xrayspectrum(energy)], [sum(dc)])
        spectrum4 = self._spectrum([c3.xrayspectrum(energy)], [sum(dc)])

        self._spectrum_equal(spectrum3, spectrum4)
        self._spectrum_equal(spectrum, spectrum3)

    @unittest.skipIf(
        xrfdetectors.compoundfromname.xraylib is None, "xraylib not installed"
    )
    def test_refractiveindex(self):
        c1 = compoundfromformula.CompoundFromFormula("Co2O3", 1.5)
        c2 = compoundfromformula.CompoundFromFormula("Fe2O3", 1.6)
        c = mixture.Mixture([c1, c2], [2, 3], types.fraction.mole)
        energy = 30

        wavelength = ureg.Quantity(energy, "keV").to("cm", "spectroscopy")
        beta = c.refractive_index_beta(energy)
        m = (
            (4 * np.pi / (wavelength * ureg.Quantity(c.density, "g/cm^3")))
            .to("cm^2/g")
            .magnitude
        )
        np.testing.assert_allclose(beta * m, c.mass_abs_coeff(energy), rtol=1e-2)

        delta = c.refractive_index_delta(energy)
        m = (
            (
                2
                * np.pi
                / (
                    ureg.classical_electron_radius
                    * wavelength ** 2
                    * ureg.particles_per_mol
                    / ureg.Quantity(c.molarmasseff(), "g/mol")
                    * c.Zeff
                )
            )
            .to("g/cm^3")
            .magnitude
        )
        np.testing.assert_allclose(delta * m, c.density, rtol=1e-2)

    @unittest.skipIf(
        xrfdetectors.compoundfromname.xraylib is None, "xraylib not installed"
    )
    def test_serialize(self):
        c1 = compoundfromformula.CompoundFromFormula("Co2O3", 1.5)
        c2 = compoundfromformula.CompoundFromFormula("Fe2O3", 1.6)
        m1 = mixture.Mixture([c1, c2], [2, 3], types.fraction.mole, name="mixture")
        m2 = jsonpickle.loads(jsonpickle.dumps(m1))
        self.assertEqual(set(m1.parts), set(m2.parts))
        self.assertEqual(m1.density, m2.density)


def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_mixture("test_molefractions"))
    testSuite.addTest(test_mixture("test_addcompounds"))
    testSuite.addTest(test_mixture("test_tocompound"))
    testSuite.addTest(test_mixture("test_cross_sections"))
    testSuite.addTest(test_mixture("test_refractiveindex"))
    testSuite.addTest(test_mixture("test_serialize"))
    return testSuite


if __name__ == "__main__":
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
