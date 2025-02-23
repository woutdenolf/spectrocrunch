import unittest
import numpy as np
import os

from testfixtures import TempDirectory
import matplotlib.pyplot as plt

from .. import pymca
from .. import multilayer
from .. import xrayspectrum
from ...geometries import xrf as xrfgeometries
from ...sources import xray as xraysources
from ...detectors import xrf as xrfdetectors
from . import xrf_setup


class test_pymca(unittest.TestCase):
    def setUp(self):
        self.dir = TempDirectory()
        self.energy = [7.5, 8]

    def tearDown(self):
        self.dir.cleanup()

    @unittest.skipIf(
        xrfdetectors.compoundfromname.xraylib is None, "xraylib not installed"
    )
    def test_loadcfg(self):
        cfgfile = os.path.join(self.dir.path, "mca.cfg")

        h1 = xrf_setup.simple(energy=self.energy)
        h1.savepymca(cfgfile)

        source = xraysources.factory("synchrotron")
        detector = xrfdetectors.factory("XRFDetector", ehole=3.8)
        geometry = xrfgeometries.factory(
            "LinearXRFGeometry",
            detector=detector,
            source=source,
            zerodistance=0,
            detectorposition=0,
            positionunits="mm",
        )
        sample = multilayer.Multilayer(geometry=geometry)

        h2 = pymca.PymcaHandle(sample=sample)
        h2.loadfrompymca(cfgfile)
        np.testing.assert_allclose(h1.mca(), h2.mca())

    @unittest.skipIf(
        xrfdetectors.compoundfromname.xraylib is None, "xraylib not installed"
    )
    def test_rates(self):
        h = xrf_setup.simple(energy=self.energy, escape=0, snip=False, continuum=1)

        if False:
            path = "/data/id21/inhouse/wout/tmp/pymcatst"
            cfgfile = os.path.join(path, "mca_mixture.cfg")
            mcafile = os.path.join(path, "spectrum.mca")
            h.savepymca(cfgfile)
            h.savemca(mcafile, func=lambda x: x + 1)

        y = h.mca()
        y += 1

        # Fit data
        h.setdata(y)
        h.configurepymca()
        fitresult = h.fit()

        # Mass fractions are calculated as follows:
        #   grouparea = flux.time.grouprate
        #   grouprate = solidangle/(4.pi).sum_l[massfrac_l.grouprate_l]  where l loops over the layers
        #
        #   massfrac_l = self.mcafit._fluoRates[layer][element]["rates"][group]      where layer in 1,...,n
        #   grouprate_l = self.mcafit._fluoRates[layer][element]["mass fraction"]
        #   sum_l[massfrac_l.grouprate_l] = self.mcafit._fluoRates[0][element]["rates"][group]*
        #                                   self.mcafit._fluoRates[0][element]["mass fraction"]
        #
        # When element in one layer:
        #   massfrac = area/(flux.time.solidangle/(4.pi).grouprate_l)
        #
        # When element in more than one layer:
        #   grouprate_avg = solidangle/(4.pi).massfrac_avg.sum_l[grouprate_l]
        #
        # When element in more than one layer (per layer as if all intensity came from that layer?):
        #   massfrac_l = grouparea/(flux.time.solidangle/(4.pi).grouprate_l)

        grouprates = h.xraygrouprates(scattering=False, method="fisx")
        safrac = h.sample.geometry.solidangle / (4 * np.pi)
        np.testing.assert_allclose(safrac, h._pymcainternals_solidanglefrac())

        for group in fitresult["fitareas"]:
            if not isinstance(group, xrayspectrum.FluoZLine):
                continue
            element, linegroup = group.element, group.group

            grouprate_avg = h.mcafit._fluoRates[0][element]["rates"][
                "{} xrays".format(linegroup)
            ]
            grouprate_avg *= safrac
            massfrac_avg = h.mcafit._fluoRates[0][element]["mass fraction"]

            grouprate = 0.0
            grouprate_avg2 = 0.0
            massfrac_avg2 = 0.0
            npresent = 0
            for j in range(h.sample.nlayers):
                i = j + 1
                if element in h.mcafit._fluoRates[i]:
                    massfrac_l = h.mcafit._fluoRates[i][element]["mass fraction"]
                    grouprate_l = h.mcafit._fluoRates[i][element]["rates"][
                        "{} xrays".format(linegroup)
                    ]
                    grouprate += massfrac_l * grouprate_l
                    grouprate_avg2 += grouprate_l
                    # massfrac_avg2 = max(massfrac_avg2,massfrac_l)
                    massfrac_avg2 = massfrac_l  # just the last one?
                    npresent += 1
            grouprate *= safrac
            grouprate_avg2 *= safrac
            # TODO: something wrong with mixtures of different elements (fisx vs Elements)
            if group in grouprates:
                # 1% error fisx vs. Elements?
                np.testing.assert_allclose(grouprate, grouprates[group], rtol=1e-2)

            np.testing.assert_allclose(grouprate, fitresult["rates"][group])

            if npresent == 1:
                np.testing.assert_allclose(grouprate_avg * massfrac_avg, grouprate)
            else:
                np.testing.assert_allclose(massfrac_avg, massfrac_avg2)
                np.testing.assert_allclose(grouprate_avg, grouprate_avg2)

            np.testing.assert_allclose(
                fitresult["massfractions"][group],
                fitresult["fitareas"][group] / (h.I0 * grouprate_avg),
            )

            for j in range(h.sample.nlayers):
                i = j + 1
                if element in h.mcafit._fluoRates[i]:
                    grouprate_l = h.mcafit._fluoRates[i][element]["rates"][
                        "{} xrays".format(linegroup)
                    ]
                    grouprate = grouprate_l * safrac
                    np.testing.assert_allclose(
                        fitresult["lmassfractions"][j][group],
                        fitresult["fitareas"][group] / (h.I0 * grouprate),
                    )

        # Plot
        plt.plot(fitresult["energy"], fitresult["y"], label="data")
        plt.plot(fitresult["energy"], fitresult["yfit"], label="pymca")

        spectrum = h.xrayspectrum()
        spectrum.plot(
            fluxtime=h.I0,
            histogram=True,
            ylog=False,
            decompose=False,
            backfunc=lambda x: 1,
        )

        ax = plt.gca()
        ax.set_ylim(ymin=np.nanmin(y[np.nonzero(y)]))
        ax.set_xlabel("Energy (keV)")
        ax.set_ylabel("Intensity (cts)")
        plt.legend(loc="best")
        # plt.show()

    @unittest.skipIf(
        xrfdetectors.compoundfromname.xraylib is None, "xraylib not installed"
    )
    def test_spectrum(self):
        h = xrf_setup.complex(
            energy=self.energy,
            escape=0,
            flux=1e10,
            time=1,
            scatter=np.zeros_like(self.energy),
            linear=1,
            emin=2,
            emax=7.4,
        )
        h.sample.geometry.detector.bltail = False
        h.sample.geometry.detector.bstail = False
        h.sample.geometry.detector.bstep = False
        y = h.mca()

        # Prepare fit
        h.setdata(y)
        h.configurepymca()

        # Force config
        config = h.mcafit.getConfiguration()
        config["fit"]["stripflag"] = 0
        h.mcafit.configure(config)

        # Fit
        # h.fitgui(loadfromfit=False)
        fitresult = h.fit(loadfromfit=False)

        # Get fundamental MCA and fitted MCA
        spectrum = h.xrayspectrum(method="fisx", scattering=False)
        x, ysum, ylabel = spectrum.sumspectrum(fluxtime=h.I0, histogram=True)
        ypymca = fitresult["interpol_energy"](fitresult["ymatrix"])(x)

        # TODO: Doesn't work due to peak rejection, add it the xrayspectrum
        # np.testing.assert_allclose(ysum,ypymca,rtol=1e-2)

        # Plot
        # x,ygroup,ylabel,names = spectrum.linespectra(fluxtime=h.I0,histogram=True)
        # for name,y in zip(names,ygroup.T):
        #    if "Ce-L3" not in name: # pymca is cutting small peaks
        #        continue
        #    plt.plot(x,y,label=name)

        plt.plot(x, ysum, label="fisx", linewidth=2)
        plt.plot(x, ypymca, label="pymca", linewidth=2)
        plt.legend()
        ax = plt.gca()
        ax.set_yscale("log", basey=10)
        plt.ylim([0.001, max(ysum)])
        # plt.show()
