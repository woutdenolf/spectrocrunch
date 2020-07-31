# -*- coding: utf-8 -*-

import unittest
import numpy as np
from scipy.spatial.transform import Rotation
from testfixtures import TempDirectory
from .. import xrmc
from .. import xrmc_utils
from ...materials import element

HAS_XRMC = xrmc.installed()


def cartesianToSpherical(p):
    x, y, z = p
    r = np.linalg.norm(p)
    if r == 0:
        polar, azimuth = 0, 0
    else:
        polar = np.degrees(np.arccos(z / r))
        azimuth = np.degrees(np.arctan2(y, x))
    return r, polar, azimuth


def sphericalToCartesian(r, polar, azimuth):
    angle = np.radians(azimuth)
    cosa = np.cos(angle)
    sina = np.sin(angle)
    angle = np.radians(polar)
    cosp = np.cos(angle)
    sinp = np.sin(angle)
    return np.array([r * sinp * cosa, r * sinp * sina, r * cosp])


class test_xrmc(unittest.TestCase):
    def setUp(self):
        self.dir = TempDirectory()

    def tearDown(self):
        print(self.dir.path)
        # self.dir.cleanup()

    @unittest.skipIf(not HAS_XRMC, "xrmc not installed")
    def test_xrf(self):
        polar_start, polar_end, polar_nsteps = 0, 180, 30
        polar = np.linspace(polar_start, polar_end, polar_nsteps + 1)
        azimuth = 0, 45, 90
        fluoline = element.Element("Ca").fluolines("KL3")

        world = xrmc.XrmcWorldBuilder(self.dir.path)
        source, pymcahandle = xrmc_utils.define_source()
        geometry = xrmc_utils.define_xrfgeometry(source=source)
        pymcahandle.sample = xrmc_utils.define_sample(geometry, name="simple")
        xrmc_utils.world_addsource(world, pymcahandle)
        xrmc_utils.world_addsdd(
            world,
            pymcahandle,
            multiplicity=10,
            polar=polar_start,
            azimuth=azimuth[0],
            convoluted=True,
        )
        xrmc_utils.world_addsample(world, pymcahandle.sample)

        polar_step = float(polar_end - polar_start) / polar_nsteps
        result = {}
        for az in azimuth:
            world.main.removeloops()
            world.detector.azimuth = az
            world.detector.add_polarrotationloop(polar_step, polar_nsteps)
            data, info = xrmc_utils.run(
                world, (1, 10000), simulate=True, ylog=True, plot=False
            )

            result[az] = data

    def assertAlmostModEqual(self, a, b, modulo=1, **kwargs):
        n = (b - a) / float(modulo)
        if "err_msg" not in kwargs:
            kwargs["err_msg"] = "{} not equal to {} (diff = {}*{})".format(
                a, b, n, modulo
            )
        idiff = n - int(np.round(n))
        np.testing.assert_allclose(idiff, 0, **kwargs)

    def test_sphericalrotation(self):
        # indirect test of xrmc_positional_device.xrmcSphericalRotationInput
        for _ in range(500):
            p1 = np.random.uniform(-1, 1, 3)
            r, polar, azimuth = cartesianToSpherical(p1)
            Rz = Rotation.from_euler("z", azimuth, degrees=True)
            Ry = Rotation.from_euler("y", polar, degrees=True)
            p2 = (Rz * Ry).apply(np.array([0, 0, r]))
            np.testing.assert_array_almost_equal(p1, p2)

    def test_polarrotation(self):
        # indirect test of xrmc_positional_device.add_polarrotationloop
        for _ in range(500):
            p1 = np.random.uniform(-1, 1, 3)
            r1, polar1, azimuth1 = cartesianToSpherical(p1)
            angle = np.radians(azimuth1 + 90)
            u = np.cos(angle)
            v = np.sin(angle)
            polarstep = np.random.uniform(-200, 200)
            rotvecs = np.radians(polarstep) * np.array([u, v, 0])
            R = Rotation.from_rotvec(rotvecs)
            p2 = R.apply(p1).T
            r2, polar2, azimuth2 = cartesianToSpherical(p2)
            polar1 = (polar1 + polarstep) % 360
            if polar1 > 180:
                polar1 = 360 - polar1
            if polar1 < 0:
                polar1 = abs(polar1)
            self.assertAlmostEqual(r1, r2)
            self.assertAlmostModEqual(azimuth1, azimuth2, modulo=180, atol=1e-10)
            self.assertAlmostModEqual(polar1, polar2, modulo=360, atol=1e-10)

    def test_azimuthalrotation(self):
        # indirect test of xrmc_positional_device.add_azimuthalrotationloop
        for _ in range(500):
            p1 = np.random.uniform(-1, 1, 3)
            r1, polar1, azimuth1 = cartesianToSpherical(p1)
            azstep = np.random.uniform(-200, 200)
            rotvecs = np.radians(azstep) * np.array([0, 0, 1])
            R = Rotation.from_rotvec(rotvecs)
            p2 = R.apply(p1).T
            r2, polar2, azimuth2 = cartesianToSpherical(p2)
            azimuth1 += azstep
            self.assertAlmostEqual(r1, r2)
            self.assertAlmostModEqual(azimuth1, azimuth2, modulo=360, atol=1e-10)
            self.assertAlmostEqual(polar1, polar2)


def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_xrmc("test_xrf"))
    # testSuite.addTest(test_xrmc("test_sphericalrotation"))
    # testSuite.addTest(test_xrmc("test_polarrotation"))
    # testSuite.addTest(test_xrmc("test_azimuthalrotation"))
    return testSuite


if __name__ == "__main__":
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
