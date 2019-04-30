import os
import numpy as np
from ...geometries import qxrf
from ...materials import compoundfromname
from ...materials import element
from ...materials import mixture
from ...materials import types
from ...materials import multilayer
from ...materials import pymca
from ...utils import listtools
from ...io import xiaedf


class XrfMapGenerator(object):

    def __init__(self, nmaps=2):
        water = compoundfromname.compoundfromname("water")
        calcite = compoundfromname.compoundfromname("calcite")
        fe = element.Element("Fe")
        mix = mixture.Mixture([water, calcite, fe], [0.5, 0.2, 0.3],
                              types.fraction.mass, name="sample")
        sample1 = multilayer.Multilayer(mix, 1e-4, geometry=None)
        mix = mixture.Mixture([water, calcite, fe], [0.3, 0.3, 0.4],
                              types.fraction.mass, name="sample")
        sample2 = multilayer.Multilayer(mix, 1e-4, geometry=None)
        pymcahandle = pymca.PymcaHandle(energy=8.0, flux=1e9, time=0.3,
                                        snip=False, continuum=1, escape=False,
                                        linear=True, noisepropagation=False)
        self.instrument = 'sxm'
        self.sample1 = sample1
        self.sample2 = sample2
        self.pymcahandle = pymcahandle
        self.energies = np.linspace(8, 9, nmaps)
        self.nmaps = nmaps
        self.ndet = 3
        self.detcomb = (0,), (1,), (2,), (0, 2), (0, 1, 2)

    def taskparams_pymca(self, seldetectors):
        parameters = self._quantparams.copy()
        parameters['xrf_positions'] = [self.xrfspectra[k]['detectorposition']
                                       for k in seldetectors]
        return parameters

    def taskparams_geometry(self, seldetectors):
        return self._taskparams_geometry(ndet=len(seldetectors))

    def _taskparams_geometry(self, ndet=1):
        init = {'reference': self.pymcahandle.flux,
                'defaultexpotime': self.pymcahandle.time,
                'xrfgeometries': [('sxm120', 'leia')]*ndet,
                'instrument': self.instrument}
        fixed = {'plot': False,
                 'dark': False,
                 'time': 1,
                 'gaindiodeI0': 1e8,
                 'gaindiodeIt': 1e7}
        variable = [{'I0_counts': 300, 'It_counts': 30,
                     'dark': True},
                    {'I0_counts': 400000, 'It_counts': 100000,
                     'energy': min(self.energies)},
                    {'I0_counts': 200000, 'It_counts': 100000,
                     'energy': max(self.energies)}]
        parameters = {'geometry': 'sxm1',
                      'fixed': fixed,
                      'init': init,
                      'variable': variable}
        return parameters

    def _instantiate_geometry(self):
        parameters = self._taskparams_geometry()
        geometry = qxrf.factory(parameters['geometry'], **parameters['init'])
        geometry.batchcalibrate_diodes(parameters['fixed'], parameters['variable'])
        geometry.diodeI0.gain *= 10
        geometry.diodeIt.gain *= 10

        self._quantparams = quantparams = {}
        quantparams['referenceflux'] = geometry.reference
        quantparams['referencetime'] = geometry.defaultexpotime
        quantparams['diodeI0gain'] = geometry.diodeI0.gain
        quantparams['diodeItgain'] = geometry.diodeI0.gain
        return geometry

    def _generate_spectra(self, path, radix):
        pymcahandle = self.pymcahandle
        samult = np.linspace(1, 0.5, self.ndet)
        bkgs = [1]*self.ndet

        xrfspectra = {}
        geometry = self._instantiate_geometry()
        geometry = geometry.xrfgeometries[0]
        solidangle = geometry.solidangle*0.6

        for seldetectors in self.detcomb:
            geometry.solidangle = solidangle*sum(samult[i] for i in seldetectors)
            bkg = sum(bkgs[i] for i in seldetectors)

            cfgfile = os.path.join(path, 'detector{}.cfg'.format(
                '_'.join(map(str, seldetectors))))
            data = {'cfgfile': cfgfile,
                    'detectorposition': geometry.detectorposition,
                    'mca': [], 'peakareas': [], 'massfractions': []}
            xrfspectra[seldetectors] = data

            for en in self.energies:
                pymcahandle.set_source(en)
                pymcahandle.emax = en+1
                mcas = []
                peakareas = []
                massfractions = []
                data['mca'].append(mcas)
                data['peakareas'].append(peakareas)
                data['massfractions'].append(massfractions)

                for sample in [self.sample1, self.sample2]:
                    sample.geometry = geometry
                    pymcahandle.sample = sample
                    pymcahandle.addtopymca(fresh=True)
                    pymcahandle.savepymca(cfgfile)

                    # pymcahandle.loadfrompymca(cfgfile)
                    #pymcahandle.emax = en+1
                    # pymcahandle.savepymca(cfgfile+'_')
                    #from ...materials.pymcadiff import diff
                    # diff(cfgfile,cfgfile+'_')
                    # exit()

                    mca = pymcahandle.mca()+bkg
                    mcas.append(mca)

                for mca in mcas:
                    pymcahandle.setdata(mca)
                    fitresult = pymcahandle.fit(loadfromfit=False)
                    # fitresult["plot"]()
                    peakareas.append(fitresult["fitareas"])
                    massfractions.append(fitresult["massfractions"])

                    labels = []
                    for label in fitresult["fitareas"]:
                        if label == 'Rayleigh':
                            labels += ["Scatter-Peak{:03d}".format(i)
                                       for i in range(label.nenergy)]
                        elif label == 'Compton':
                            labels += ["Scatter-Compton{:03d}".format(i)
                                       for i in range(label.nenergy)]
                        else:
                            labels.append(str(label))

        self.include_detectors = [(0, 2), 2, None, (1, (0, 2))]
        self.xrfspectra = xrfspectra
        self.labels = labels
        self.nchan = len(mca)
        self._check_spectra()

    def _check_spectra(self):
        # check detector sum/average
        def farr(x):
            return listtools.numpy_flatten(x)
        xrfspectra = self.xrfspectra
        energies = self.energies
        for seldetectors, data in xrfspectra.items():
            for i in range(len(energies)):
                for j in range(2):
                    a = sum(xrfspectra[(k,)]['mca'][i][j]
                            for k in seldetectors)
                    b = data['mca'][i][j]
                    np.testing.assert_allclose(a, b)
                    a = sum(farr(xrfspectra[(k,)]['peakareas'][i][j].values())
                            for k in seldetectors)
                    b = farr(data['peakareas'][i][j].values())
                    np.testing.assert_allclose(a, b)
                    a = sum(farr(xrfspectra[(k,)]['massfractions'][i][j].values())
                            for k in seldetectors)
                    b = farr(data['massfractions'][i][j].values())
                    np.testing.assert_allclose(a/len(seldetectors), b)

    def generate(self, path, radix, applyflux=True, applydt=True):
        self._generate_spectra(path, radix)
        qxrfgeometry = self._instantiate_geometry()
        refflux = qxrfgeometry.reference.to("hertz").magnitude
        reftime = qxrfgeometry.defaultexpotime.to("seconds").magnitude

        # 3 maps of a moving hotspot
        nlines, nspec = 7, 6
        nmaps = self.nmaps
        ndet = self.ndet
        nchan = self.nchan
        data = np.ones((nmaps, nlines, nspec, nchan, ndet))
        off = max(min(nlines, nspec)-nmaps, 0)//2
        for idet in range(ndet):
            det = self.xrfspectra[(idet,)]
            for imap, spectra in enumerate(det['mca']):
                spec1, spec2 = spectra
                data[imap, ..., idet] = spec1
                data[imap, imap+off, imap+off, :, idet] = spec2

        # Generate counter headers
        elabel = qxrfgeometry.instrument.edfheaderkeys['energylabel']
        tlabel = qxrfgeometry.instrument.edfheaderkeys['timelabel']
        ctrheaders = np.vectorize(lambda e: {elabel: e, tlabel: reftime},
                                  otypes=[object])(self.energies)

        # Init counters
        ctrs = {}

        # Apply flux decay
        if applyflux:
            rflux = np.linspace(1, 0.5, nmaps*nlines *
                                nspec).reshape((nmaps, nlines, nspec))
        else:
            rflux = np.ones((nmaps, nlines, nspec), dtype=float)
        flux = rflux*refflux
        fluxT = flux  # TODO real transmission
        data *= rflux[..., np.newaxis, np.newaxis]

        ctrs["arr_iodet"] = np.ones((nmaps, nlines, nspec))
        ctrs["arr_idet"] = np.ones((nmaps, nlines, nspec))
        ctrs["arr_norm"] = np.ones((nmaps, nlines, nspec))
        for i, en in enumerate(self.energies):
            ctrs["arr_iodet"][i, ...] = qxrfgeometry.diodeI0.fluxtocps(
                en, flux[i, ...])*reftime
            ctrs["arr_idet"][i, ...] = qxrfgeometry.diodeIt.fluxtocps(
                en, fluxT[i, ...])*reftime
            ctrs["arr_norm"][i, ...] = rflux[i, ...]
            op, fref, tref, traw = qxrfgeometry.xrfnormop(en)
            assert(fref == refflux)
            assert(tref == reftime)
            np.testing.assert_allclose(
                ctrs["arr_norm"][i, ...], op(ctrs["arr_iodet"][i, ...]))
            op, _ = qxrfgeometry.I0op(en, expotime=reftime)
            np.testing.assert_allclose(
                flux[i, ...], op(ctrs["arr_iodet"][i, ...]))
            op, _ = qxrfgeometry.Itop(en, expotime=reftime)
            np.testing.assert_allclose(
                flux[i, ...], op(ctrs["arr_idet"][i, ...]))

        # Deadtime corrected counters
        a = nchan//2
        b = a+100
        for i in range(ndet):
            ctrs["xmap_x1c_{:02d}".format(i)] = data[..., a:b, i].sum(axis=-1)

        b = a
        a = a-100
        for i in range(ndet):
            ctrs["xmap_x2c_{:02d}".format(i)] = data[..., a:b, i].sum(axis=-1)

        # Apply deadtime
        stats = np.zeros(
            (nmaps, nlines, nspec, xiaedf.xiadata.NSTATS, ndet), dtype=data.dtype)
        RT = reftime  # sec
        for i in range(ndet):
            Rreal = data[..., i].sum(axis=-1)
            if applydt:
                DTslow = 0.2*i/(ndet-1.)
            else:
                DTslow = 0
            DTfast = 0.05*DTslow
            LTfast = RT*(1-DTfast)
            Rslow = Rreal*(1-DTslow)
            Nslow = RT*Rslow
            ctrs["xmap_icr_{:02d}".format(i)] = Rreal
            ctrs["xmap_ocr_{:02d}".format(i)] = Rslow
            data[..., i] = data[..., i]*(1-DTslow)[..., np.newaxis]
            stats[..., xiaedf.xiadata.STDET, i] = i
            stats[..., xiaedf.xiadata.STEVT, i] = Nslow
            stats[..., xiaedf.xiadata.STICR, i] = Rreal
            stats[..., xiaedf.xiadata.STOCR, i] = Rslow
            stats[..., xiaedf.xiadata.STDT, i] = DTslow*100  # %
            stats[..., xiaedf.xiadata.STLT, i] = LTfast*1000  # msec

        # Generate data
        stack = xiaedf.xiastack_radix(path, radix)
        xialabels = ["xia{:02d}".format(i) for i in range(ndet)]
        stack.save(data, xialabels, stats=stats,
                   ctrs=np.stack(tuple(ctrs.values()), axis=-1),
                   ctrnames=ctrs.keys(),
                   ctrheaders=ctrheaders)
        
        self.path = path
        self.radix = radix
        self.data = data
        self.stats = stats
        self.ctrs = ctrs

    @property
    def scannumbers(self):
        return list(range(self.nmaps))
