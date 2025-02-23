import unittest

from .. import compound as compoundraw
from .. import compoundfromformula
from .. import compoundfromlist
from .. import compoundfromcif
from .. import compoundfromname
from .. import types
from ...patch.pint import ureg
from ...patch import jsonpickle

import numpy as np
import warnings

try:
    import iotbx.cif as iotbxcif
except ImportError:
    iotbxcif = None
    warnings.warn("cctbx is not installed", ImportWarning)

try:
    import xraylib
except ImportError:
    xraylib = None
    warnings.warn("xraylib is not installed", ImportWarning)


class test_compound(unittest.TestCase):
    @unittest.skipIf(xraylib is None, "xraylib not installed")
    def test_comparable(self):
        c1 = compoundfromformula.CompoundFromFormula(
            "C6H2(NO2)3CH3", 1.2, name="compound"
        )
        c2 = compoundfromformula.CompoundFromFormula(
            "C6H2(NO2)3CH3", 1.2, name="compound"
        )
        c3 = compoundfromformula.CompoundFromFormula(
            "C6H2(NO2)3CH3", 1.21, name="compound"
        )
        c4 = compoundfromformula.CompoundFromFormula(
            "C6H2(NO2)3CH3", 1.2, name="compoundx"
        )
        c5 = compoundfromformula.CompoundFromFormula(
            "C6H2(NO2)3CH2", 1.2, name="compound"
        )
        self.assertEqual(c1, c2)
        self.assertEqual(c1, c3)
        self.assertNotEqual(c1, c4)
        self.assertEqual(c1, c5)  # this is by design but may be unwanted?

    @unittest.skipIf(xraylib is None, "xraylib not installed")
    def test_formula(self):
        elements = ["C", "N", "O", "H"]
        a = [7, 3, 6, 5]
        density = 2.3
        c = compoundfromformula.CompoundFromFormula("C6H2(NO2)3CH3", density)
        elements2 = c.equivalents()
        for i in range(len(elements)):
            self.assertEqual(elements2[elements[i]], a[i])
        self.assertEqual(c.density, density)

    @unittest.skipIf(xraylib is None, "xraylib not installed")
    def test_list(self):
        elements = ["Fe", "S", "O"]
        a = [1, 1, 4.0]
        density = 2.3
        c = compoundfromlist.CompoundFromList(elements, a, types.fraction.mole, density)
        elements2 = c.equivalents()
        for i in range(len(elements)):
            self.assertEqual(elements2[elements[i]], a[i])
        self.assertEqual(c.density, density)
        c = compoundfromlist.CompoundFromList(elements, a, types.fraction.mass, density)
        wfrac = c.massfractions()
        for i in range(len(elements)):
            self.assertAlmostEqual(wfrac[elements[i]], a[i] / float(sum(a)))

    @unittest.skipIf(iotbxcif is None, "cctbx not installed")
    def test_cif(self):
        elements = ["Ca", "C", "O"]
        a = [6, 6, 18.0]  # unit cell content
        c = compoundfromcif.CompoundFromCif("cif/calcite.cif", name="calcite")
        elements2 = c.equivalents()
        for i in range(len(elements)):
            self.assertEqual(elements2[elements[i]], a[i])

    @unittest.skipIf(xraylib is None, "xraylib not installed")
    def test_addelements(self):
        c1 = compoundraw.Compound(["Fe", "O"], [2, 3], types.fraction.mole, density=1)
        n1 = c1.molefractions()
        snfrac1 = sum(c1.equivalents().values())
        for fractype in ["mass", "mole"]:
            c2 = compoundraw.Compound(["Fe"], [5], types.fraction.mole, density=1)
            if fractype == "mass":
                c2.addelements("O", c1.massfractions()["O"], types.fraction.mass)
            else:
                c2.addelements("O", n1["O"], types.fraction.mole)
            n2 = c2.molefractions()
            snfrac2 = sum(c2.equivalents().values())
            self.assertEqual(set(n1.keys()), set(n2.keys()))
            if fractype == "mole":
                self.assertEqual(snfrac1, snfrac2)
            else:
                self.assertEqual(1, snfrac2)
            for k in n1:
                self.assertAlmostEqual(n1[k], n2[k])

    def test_vacuum(self):
        c = compoundraw.Compound([], [], types.fraction.mole, 0, name="vacuum")
        self.assertEqual(len(c.elements), 0)
        self.assertEqual(len(c.massfractions()), 0)
        self.assertEqual(len(c.molefractions()), 0)
        self.assertEqual(c.molarmass(), 0)
        self.assertEqual(c.density, 0)

    @unittest.skipIf(xraylib is None, "xraylib not installed")
    def test_name(self):
        c = compoundfromname.compoundfromname("vacuum")
        self.assertEqual(len(c.elements), 0)
        self.assertEqual(len(c.massfractions()), 0)
        self.assertEqual(len(c.molefractions()), 0)
        self.assertEqual(c.molarmass(), 0)
        self.assertEqual(c.density, 0)

        c = compoundfromname.compoundfromname("linseed oil")

        with self.assertRaises(KeyError):
            c = compoundfromname.compoundfromname("linseed oill")

    @unittest.skipIf(xraylib is None, "xraylib not installed")
    def test_refractiveindex(self):
        density = float(2.328)
        c = compoundfromformula.CompoundFromFormula("Si", density, name="silicon")
        energy = np.asarray([5, 10], dtype=float)
        Z = 14

        # Xraylib (TODO: n_im sign is wrong!)
        n_re0 = np.asarray(
            [xraylib.Refractive_Index_Re("Si", e, density) for e in energy]
        )
        n_im0 = -np.asarray(
            [xraylib.Refractive_Index_Im("Si", e, density) for e in energy]
        )

        # Exactly like Xraylib (TODO: CS_Total -> CS_Photo_Total)
        delta = (
            density
            * 4.15179082788e-4
            * (Z + np.asarray([xraylib.Fi(Z, e) for e in energy]))
            / xraylib.AtomicWeight(Z)
            / energy**2
        )
        beta = (
            np.asarray([xraylib.CS_Total(Z, e) for e in energy])
            * density
            * 9.8663479e-9
            / energy
        )
        n_re1 = 1 - delta
        n_im1 = -beta

        np.testing.assert_allclose(n_re0, n_re1)
        np.testing.assert_allclose(n_im0, n_im1)

        # Im: Kissel
        n_im1b = (
            -np.asarray([xraylib.CS_Total_Kissel(Z, e) for e in energy])
            * density
            * 9.8663479e-9
            / energy
        )

        # Im: Kissel with pint
        n_im1d = (
            ureg.Quantity(c.mass_att_coeff(energy), "cm^2/g")
            * ureg.Quantity(c.density, "g/cm^3")
            * ureg.Quantity(energy, "keV").to("cm", "spectroscopy")
            / (4 * np.pi)
        )
        n_im1d = -n_im1d.to("dimensionless").magnitude

        r = n_im1b / n_im1d
        np.testing.assert_allclose(r, r[0])

        # Im: other formula
        n_im1c = (
            density
            * 4.15179082788e-4
            * (np.asarray([xraylib.Fii(Z, e) for e in energy]))
            / xraylib.AtomicWeight(Z)
            / energy**2
        )

        # Spectrocrunch
        n_re2 = c.refractive_index_real(energy)
        n_im2 = c.refractive_index_imag(energy)

        np.testing.assert_allclose(n_re2, n_re0)
        np.testing.assert_allclose(n_im2, n_im1c, rtol=1e-6)

    @unittest.skipIf(xraylib is None, "xraylib not installed")
    def test_refractiveindex2(self):
        density = 5.3
        c = compoundfromformula.CompoundFromFormula("Fe2O3", density, name="test")
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
                    * wavelength**2
                    * ureg.particles_per_mol
                    / ureg.Quantity(c.molarmasseff(), "g/mol")
                    * c.Zeff
                )
            )
            .to("g/cm^3")
            .magnitude
        )
        np.testing.assert_allclose(delta * m, c.density, rtol=1e-2)

    @unittest.skipIf(xraylib is None, "xraylib not installed")
    def test_serialize(self):
        c1 = compoundfromformula.CompoundFromFormula("Fe2O3", 5.3, name="test")
        c2 = jsonpickle.loads(jsonpickle.dumps(c1))
        self.assertEqual(c1, c2)
        self.assertEqual(set(c1.elements), set(c2.elements))
        self.assertEqual(c1.density, c2.density)
        c1 = compoundfromname.compoundfromname("ultralene")
        c2 = jsonpickle.loads(jsonpickle.dumps(c1))
        self.assertEqual(c1, c2)
        self.assertEqual(set(c1.elements), set(c2.elements))
        self.assertEqual(c1.density, c2.density)


def main_test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_compound("test_comparable"))
    testSuite.addTest(test_compound("test_formula"))
    testSuite.addTest(test_compound("test_list"))
    testSuite.addTest(test_compound("test_cif"))
    testSuite.addTest(test_compound("test_addelements"))
    testSuite.addTest(test_compound("test_vacuum"))
    testSuite.addTest(test_compound("test_name"))
    testSuite.addTest(test_compound("test_refractiveindex"))
    testSuite.addTest(test_compound("test_refractiveindex2"))
    testSuite.addTest(test_compound("test_serialize"))
    return testSuite


if __name__ == "__main__":
    import sys

    mysuite = main_test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
