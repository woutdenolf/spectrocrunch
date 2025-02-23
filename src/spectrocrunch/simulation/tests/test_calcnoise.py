import unittest

from ...math import noisepropagation
from .. import calcnoise
from ...materials import multilayer
from ...geometries import flatarea
from ...sources import xray as xraysources
from ...detectors import area
from ...materials.compoundfromformula import CompoundFromFormula as compound

import numpy as np


class test_calcnoise(unittest.TestCase):
    @unittest.skipIf(calcnoise.lenses.visirlib.PyTMM is None, "PyTMM not installed")
    @unittest.skipIf(
        multilayer.compoundfromdb.compoundfromname.xraylib is None,
        "xraylib not installed",
    )
    def test_ffnoise(self):
        flux = 1e5
        energy = np.linspace(3, 5, 100)
        tframe = 0.07
        nframe = 100
        ndark = 10

        src = xraysources.factory("synchrotron")
        detector = area.factory("PCO Edge 5.5")
        geometry = flatarea.factory("perpendicular", detector=detector, source=src)
        sample = multilayer.Multilayer(
            material=compound("CaCO3", 2.71), thickness=5e-4, geometry=geometry
        )
        kwargs = {
            "tframe_data": tframe,
            "nframe_data": nframe,
            "tframe_flat": tframe,
            "nframe_flat": nframe,
            "nframe_dark": ndark,
        }
        o = calcnoise.id21_ffsetup(sample=sample)
        signal, noise = o.xanes(flux, energy, **kwargs)
        self.assertEqual(signal.shape, (energy.size, 1))
        o.plotxanesnoise(flux, energy, **kwargs)
        # plt.show()

    @unittest.skipIf(calcnoise.lenses.visirlib.PyTMM is None, "PyTMM not installed")
    @unittest.skipIf(
        multilayer.compoundfromdb.compoundfromname.xraylib is None,
        "xraylib not installed",
    )
    def test_reverse(self):
        flux = np.linspace(1e4, 1e5, 2)
        energy = np.linspace(3, 5, 3)
        tframe = 0.07
        nframe = 100

        src = xraysources.factory("synchrotron")
        detector = area.factory("PCO Edge 5.5")
        geometry = flatarea.factory("perpendicular", detector=detector, source=src)
        sample = multilayer.Multilayer(
            material=compound("CaCO3", 2.71), thickness=5e-4, geometry=geometry
        )
        o = calcnoise.id21_ffsetup(sample=sample)
        ph1 = np.broadcast_to(flux * tframe, (energy.size, flux.size))
        Ndet2 = None
        for withnoise in [False, True]:
            if withnoise:
                N = noisepropagation.poisson(flux * tframe)
            else:
                N = flux * tframe
            Ndet = o.propagate(
                N, energy, tframe=tframe, nframe=nframe, forward=True, samplein=True
            )
            if withnoise and Ndet2 is not None:
                np.testing.assert_allclose(Ndet2, noisepropagation.E(Ndet))
            self.assertEqual(Ndet.shape, (energy.size, flux.size))
            ph2 = o.propagate(
                Ndet, energy, tframe=tframe, nframe=nframe, forward=False, samplein=True
            )
            self.assertEqual(ph2.shape, (energy.size, flux.size))
            if withnoise:
                np.testing.assert_allclose(ph1, noisepropagation.E(ph2))
                np.testing.assert_allclose(ph1, noisepropagation.VAR(ph2))
            else:
                np.testing.assert_allclose(ph1, ph2)
            Ndet2 = Ndet

    @unittest.skipIf(calcnoise.lenses.visirlib.PyTMM is None, "PyTMM not installed")
    @unittest.skipIf(
        multilayer.compoundfromdb.compoundfromname.xraylib is None,
        "xraylib not installed",
    )
    def test_totalgain(self):
        flux = np.linspace(1e4, 1e5, 10)
        energy = 5
        tframe = 0.07
        nframe = 100

        o = calcnoise.id21_ffsetup()
        ph = flux * tframe
        DU = o.propagate(ph, energy, tframe=tframe, nframe=nframe)
        m, b = o.photontoDU(energy, tframe, nframe)
        np.testing.assert_allclose(DU[0, :], ph * m + b)
