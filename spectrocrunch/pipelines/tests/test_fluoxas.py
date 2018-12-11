# -*- coding: utf-8 -*-
#
#   Copyright (C) 2017 European Synchrotron Radiation Facility, Grenoble, France
#
#   Principal author:   Wout De Nolf (wout.de_nolf@esrf.eu)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import unittest
import numpy as np
import os
import sys
import contextlib
import itertools
import logging
import re
from testfixtures import TempDirectory
from PyMca5.PyMcaIO import ConfigDict

from .. import fluoxas
from ..run import run_sequential
from ...io import xiaedf
from ...io import nxfs
from ...process import nxresult
from ...align import types
from ...utils import instance
from ...utils import listtools
from ...geometries import qxrf
from ...geometries import xrf as xrfgeometries
from ...sources import xray as xraysources
from ...detectors import xrf as xrfdetectors
from ...materials import compoundfromname
from ...materials import compoundfromformula
from ...materials import element
from ...materials import mixture
from ...materials import types
from ...materials import multilayer
from ...materials import pymca

try:
    logger = logging.getLogger(__loader__.fullname)
except NameError:
    logger = logging.getLogger(__name__)
    
class test_fluoxas(unittest.TestCase):

    def setUp(self):
        self.dir = TempDirectory()
    
    def tearDown(self):
        self.dir.cleanup()
        
    def _data_init(self):
        synchrotron = xraysources.factory("synchrotron")
        leia = xrfdetectors.factory("leia")
        geometry = None
        water = compoundfromname.compoundfromname("water")
        calcite = compoundfromname.compoundfromname("calcite")
        fe = element.Element("Fe")
        mix = mixture.Mixture([water,calcite,fe],[0.5,0.2,0.3],types.fraction.mass,name="sample")
        sample1 = multilayer.Multilayer(mix,1e-4,geometry=geometry)
        mix = mixture.Mixture([water,calcite,fe],[0.3,0.3,0.4],types.fraction.mass,name="sample")
        sample2 = multilayer.Multilayer(mix,1e-4,geometry=geometry)
        pymcahandle = pymca.PymcaHandle(energy=8.0,flux=1e9,time=0.1,snip=False,
                    continuum=1,escape=False,linear=True,noisepropagation=False)

        energy = [8.,9.]
        
        ndet = 3
        samult = [1,0.8,0.5]
        bkgs = [1,1,1]

        detectors = {}
        for seldetectors in [(0,),(1,),(2,),(0,2),(0,1,2)]:
            geometry = xrfgeometries.factory("sxm120",detector=leia,
                                             source=synchrotron,detectorposition=-15.)
            geometry.solidangle *= sum(samult[i] for i in seldetectors)
            bkg = sum(bkgs[i] for i in seldetectors)
            
            cfgfile = os.path.join(self.dir.path,'detector{}.cfg'.format('_'.join(map(str,seldetectors))))
            detector = {'cfgfile':cfgfile,'geometry':geometry,
                        'mca':[],'peakareas':[],'massfractions':[]}
            detectors[seldetectors] = detector
            
            for en in energy:
                pymcahandle.set_source(en)
                pymcahandle.emax = en+1
                mcas = []
                peakareas = []
                massfractions = []
                detector['mca'].append(mcas)
                detector['peakareas'].append(peakareas)
                detector['massfractions'].append(massfractions)
                
                for sample in [sample1,sample2]:
                    sample.geometry = geometry
                    pymcahandle.sample = sample
                    pymcahandle.addtopymca(fresh=True)
                    pymcahandle.savepymca(cfgfile)
                    
                    #pymcahandle.loadfrompymca(cfgfile)
                    #pymcahandle.emax = en+1
                    #pymcahandle.savepymca(cfgfile+'_')
                    #from ...materials.pymcadiff import diff
                    #diff(cfgfile,cfgfile+'_')
                    #exit()
                    
                    mca = pymcahandle.mca()+bkg
                    mcas.append(mca)
 
                for mca in mcas:
                    pymcahandle.setdata(mca)
                    fitresult = pymcahandle.fit(loadfromfit=False)
                    #fitresult["plot"]()
                    peakareas.append(fitresult["fitareas"])
                    massfractions.append(fitresult["massfractions"])

                    labels = []
                    for label in fitresult["fitareas"]:
                        if label=='Rayleigh':
                            labels += ["Scatter-Peak{:03d}".format(i)
                                        for i in range(label.nenergy)]
                        elif label=='Compton':
                            labels += ["Scatter-Compton{:03d}".format(i)
                                        for i in range(label.nenergy)]
                        else:
                            labels.append(str(label))
                        
        self.procinfo = {}
        self.procinfo['include_detectors'] = [(0,2),2,None,(1,(0,2))]
        self.procinfo['flux'] = pymcahandle.flux
        self.procinfo['time'] = pymcahandle.time
        self.procinfo['detectors'] = detectors
        self.procinfo['energy'] = energy
        self.procinfo['ndet'] = ndet
        self.procinfo['nchan'] = len(mca)
        self.procinfo['labels'] = labels
        
        # check detector sum/average
        for seldetectors,detector in detectors.items():
            for i in range(len(energy)):
                for j in range(2):
                    a = sum(detectors[(k,)]['mca'][i][j] for k in seldetectors)
                    b = detector['mca'][i][j]
                    np.testing.assert_allclose(a,b)
                    
                    a = sum(np.array(detectors[(k,)]['peakareas'][i][j].values())
                            for k in seldetectors)
                    b = np.array(detector['peakareas'][i][j].values())
                    np.testing.assert_allclose(a,b)
                    
                    a = sum(np.array(detectors[(k,)]['massfractions'][i][j].values())
                            for k in seldetectors)
                    b = np.array(detector['massfractions'][i][j].values())
                    np.testing.assert_allclose(a/len(seldetectors),b)

    def _data_generate(self,applyflux=True,applydt=True):
        self._data_init()
        qxrfgeometry = self.qxrfgeometry()
        refflux = qxrfgeometry.reference.to("hertz").magnitude
        expotime = qxrfgeometry.defaultexpotime.to("seconds").magnitude

        # 3 maps of a moving hotspot
        nlines,nspec = 7,6
        nmaps = len(self.procinfo['energy'])
        ndet = self.procinfo['ndet']
        nchan = self.procinfo['nchan']
        data = np.ones((nmaps,nlines,nspec,nchan,ndet))
        off = max(min(nlines,nspec)-nmaps,0)//2
        for idet in range(ndet):
            det = self.procinfo['detectors'][(idet,)]
            for imap,spectra in enumerate(det['mca']):
                spec1,spec2 = spectra
                data[imap,...,idet] = spec1
                data[imap,imap+off,imap+off,:,idet] = spec2

        # Generate counter headers
        energy = self.procinfo['energy']
        ctrheaders = np.vectorize(lambda e:{"DCM_Energy":e,"time":expotime},
                                  otypes=[object])(energy)
        
        # Init counters
        ctrs = {}
        
        # Apply flux decay
        if applyflux:
            rflux = np.linspace(1,0.5,nmaps*nlines*nspec).reshape((nmaps,nlines,nspec))
        else:
            rflux = np.ones((nmaps,nlines,nspec),dtype=float)
        flux = rflux*refflux
        data *= rflux[...,np.newaxis,np.newaxis]

        ctrs["arr_iodet"] = np.ones((nmaps,nlines,nspec))
        ctrs["arr_idet"] = np.ones((nmaps,nlines,nspec))
        ctrs["arr_norm"] = np.ones((nmaps,nlines,nspec))
        for i,en in enumerate(energy):
            ctrs["arr_iodet"][i,...] = qxrfgeometry.diodeI0.fluxtocps(en,flux[i,...])*expotime
            ctrs["arr_idet"][i,...] = qxrfgeometry.diodeIt.fluxtocps(en,flux[i,...])*expotime #TODO real transmission
            ctrs["arr_norm"][i,...] = rflux[i,...]
            op,fref,tref,traw = qxrfgeometry.xrfnormop(en)
            self.assertEqual(fref,refflux)
            np.testing.assert_allclose(ctrs["arr_norm"][i,...],op(ctrs["arr_iodet"][i,...]))
            op,_ = qxrfgeometry.I0op(en,expotime=expotime)
            np.testing.assert_allclose(flux[i,...],op(ctrs["arr_iodet"][i,...]))
            op,_ = qxrfgeometry.Itop(en,expotime=expotime)
            np.testing.assert_allclose(flux[i,...],op(ctrs["arr_idet"][i,...]))

        # Deadtime corrected counters
        a = nchan/2
        b = a+100
        for i in range(ndet):
            ctrs["xmap_x1c_{:02d}".format(i)] = data[...,a:b,i].sum(axis=-1)

        b = a
        a = a-100
        for i in range(ndet):
            ctrs["xmap_x2c_{:02d}".format(i)] = data[...,a:b,i].sum(axis=-1)

        # Apply deadtime
        stats = np.zeros((nmaps,nlines,nspec,xiaedf.xiadata.NSTATS,ndet),dtype=data.dtype)
        for i in range(ndet):
            ICR = data[...,i].sum(axis=-1)
            if applydt:
                OCR = ICR*(1-0.2*i/(ndet-1.))# DT = 10*i %
            else:
                OCR = ICR
            ctrs["xmap_icr_{:02d}".format(i)] = ICR
            ctrs["xmap_ocr_{:02d}".format(i)] = OCR
            
            data[...,i] = data[...,i]*(OCR/ICR)[...,np.newaxis]

            stats[...,xiaedf.xiadata.STDET,i] = i
            stats[...,xiaedf.xiadata.STEVT,i] = ICR # % Not sure
            stats[...,xiaedf.xiadata.STICR,i] = ICR
            stats[...,xiaedf.xiadata.STOCR,i] = OCR
            stats[...,xiaedf.xiadata.STDT,i] = 100-OCR*100./ICR # %
            stats[...,xiaedf.xiadata.STLT,i] = OCR*1000./ICR # 1000 msec RT

        # Generate data
        path = self.dir.path
        radix = "test"
        stack = xiaedf.xiastack_radix(path,radix)
        xialabels = ["xia{:02d}".format(i) for i in range(ndet)]
        stack.save(data,xialabels,stats=stats,ctrs=np.stack(ctrs.values(),axis=-1),ctrnames=ctrs.keys(),ctrheaders=ctrheaders)

        self._data = path,radix,data,stats,ctrs,qxrfgeometry
                              
    def test_process(self):
        self._data_generate()
        
        parameters = [(None,),(True,False),self.procinfo['include_detectors'],(True,False),
                      (False,True),(True,False),(True,False),(0,),(False,True)]
        
        if hasattr(self,'subTest'):
            for i,combination in enumerate(itertools.product(*parameters)):
                with self.subTest(i=i):
                    self._process(*combination)
        else:
            for i,combination in enumerate(itertools.product(*parameters)):
                self._process(*combination)
                sys.stdout.write('.')
                sys.stdout.flush()
                    
    def _process(self,alignmethod,cfgfileuse,include_detectors_p,adddetectors_p,\
                      addbeforefit,quant,dtcor,stackdim,correctspectra):

        if not cfgfileuse and alignmethod is not None:
            return
        
        sourcepath,radix,data,stats,ctrs,qxrfgeometry = self._data
        nmaps,nlines,nspec,nchan,ndet = data.shape
        scannumbers = [range(nmaps)]

        if include_detectors_p:
            include_detectors = include_detectors_p
        else:
            include_detectors = tuple(range(ndet))
            
        alldetectors = tuple(sorted(list(listtools.flatten(include_detectors))))
        adddetectorgroups = any(len(instance.asarray(dets))>1 for dets in instance.asarray(include_detectors))
        adddetectors = adddetectors_p and len(alldetectors)>1
        addspectra = (adddetectors or adddetectorgroups) and addbeforefit
        fluxnormbefore = quant and correctspectra
        dtcorbefore = dtcor and (correctspectra or addspectra)
        newspectra = addspectra or fluxnormbefore or dtcorbefore
        
        if addspectra:
            if adddetectorgroups:
                seldetectors = [tuple(instance.asarray(dets).tolist()) for dets in instance.asarray(include_detectors)]
            else:
                seldetectors = [alldetectors]
        else:
            seldetectors = [(det,) for det in alldetectors]

        if quant:
            geom = qxrfgeometry
            geom.xrfgeometries = [self.procinfo['detectors'][k]['geometry'] for k in seldetectors]
            prealignnormcounter = None
        else:
            geom = None
            prealignnormcounter = "arr_norm"

        if cfgfileuse:
            cfgfiles = [self.procinfo['detectors'][k]['cfgfile'] for k in seldetectors]
        else:
            cfgfiles = None

        alignreference = None
        fitlabels = set(self.fitlabels(quant=quant))
        fitlabelsfile = set([label.replace('-','_') for label in fitlabels])
        detcounterlabels = set(['xmap_icr','xmap_ocr','xmap_x1c','xmap_x2c'])
        counterlabels = set(['arr_iodet','arr_idet','arr_norm'])
        calclabels = set(['calc_transmission','calc_absorbance','calc_flux0','calc_fluxt'])
        for label in fitlabels:
            if not 'Scatter' in label:
                alignreference = label
                break
        refimageindex = 0

        # Data
        expectedgroups_data = []
        if addspectra:
            if adddetectorgroups:
                expectedgroups_data = set(["S{:d}".format(i+1) for i in range(len(include_detectors))])
            else:
                if len(alldetectors)==1:
                    expectedgroups_data = set(['{:02d}'.format(alldetectors[0])])
                else:
                    expectedgroups_data = set(['S1'])
        else:
            expectedgroups_data = set(["{:02d}".format(i) for i in list(listtools.flatten(include_detectors))])
            
        # Final groups
        if adddetectors:
            if adddetectorgroups:
                expectedgroups_result = ["S{:d}".format(len(include_detectors)+1)]
            else:
                expectedgroups_result = ['S1']
        elif adddetectorgroups:
            expectedgroups_result = ["S{:d}".format(i+1) for i in range(len(include_detectors))]
        else:
            expectedgroups_result = ["{:02d}".format(i) for i in alldetectors]
        
        alignref = "/detector{}/{}".format(expectedgroups_result[0],alignreference)
        expectedgroups_result = ['counters'] + ["detector"+det for det in expectedgroups_result]
        expectedgroups_result = set(expectedgroups_result)
        
        # Processes
        expected_nxprocess = ['process:pymca.1']
        if prealignnormcounter is not None:
            expected_nxprocess.append('process:normalize.1')
        if alignmethod is not None and alignreference is not None:
            expected_nxprocess.append('process:align.1')
            expected_nxprocess.append('process:crop.1')
        
        with self._destpath_context() as destpath:
            for repeat in range(2):
                logger.debug("alignmethod = {}".format(alignmethod))
                logger.debug("stackdim = {}".format(stackdim))
                logger.debug("cfgfiles = {}".format(cfgfiles))
                logger.debug("dtcor = {}".format(dtcor))
                logger.debug("adddetectors = {}".format(adddetectors_p))
                logger.debug("addbeforefit = {}".format(addbeforefit))
                logger.debug("include_detectors = {}".format(include_detectors_p))
                logger.debug("adddetectorgroups = {}".format(adddetectorgroups))
                logger.debug("quant = {}".format(quant))
                logger.debug("correctspectra = {}".format(correctspectra))
                logger.debug("newspectra = {}".format(newspectra))
                
                parameters = {}
                parameters["nxentry"] = os.path.join(destpath.path,radix+'.h5')+':/'+radix
                parameters["sourcepath"] = sourcepath
                parameters["scanname"] = radix
                parameters["scannumbers"] = scannumbers
                parameters["cfgfiles"] = cfgfiles
                parameters["alignmethod"] = alignmethod
                parameters["alignreference"] = alignref
                parameters["refimageindex"] = refimageindex
                parameters["dtcor"] = dtcor
                parameters["plot"] = False
                parameters["adddetectors"] = adddetectors_p
                parameters["addbeforefit"] = addbeforefit
                parameters["qxrfgeometry"] = geom
                parameters["correctspectra"] = correctspectra
                parameters["replacenan"] = bool(alignmethod)
                parameters["prealignnormcounter"] = prealignnormcounter
                parameters["stackdim"] = stackdim
                parameters["include_detectors"] = include_detectors_p
                parameters["instrument"] = "id21"
                parameters["counters"] = ["arr_norm"]
                parameters["fastfitting"] = True
                parameters["replacenan"] = False
                parameters["crop"] = alignmethod is not None

                tasks = fluoxas.tasks(**parameters)
                if repeat:
                    for task in tasks:
                        self.assertTrue(task.done)
                    continue
                else:
                    for task in tasks:
                        self.assertFalse(task.done)
                    run_sequential(tasks)
                    for task in tasks:
                        self.assertTrue(task.done)
                    nxprocess = tasks[-1].output

                # Check generated spectra (files)
                if newspectra:
                    corlabel = ""
                    if dtcorbefore:
                        corlabel = corlabel+"dt"
                    if fluxnormbefore:
                        corlabel = corlabel+"fl"
                    if corlabel:
                        radixout = "{}_{}cor".format(radix,corlabel)
                    else:
                        radixout = radix

                    if addspectra:
                        expected = ["{}_xia{}_{:04d}_0000_{:04d}.edf".format(radixout,det,mapnum,linenum)\
                                                                    for det in expectedgroups_data\
                                                                    for mapnum in range(nmaps)\
                                                                    for linenum in range(nlines)]
                    else:
                        expected = ["{}_xia{}_{:04d}_0000_{:04d}.edf".format(radixout,det,mapnum,linenum)\
                                                                    for det in expectedgroups_data\
                                                                    for mapnum in range(nmaps)\
                                                                    for linenum in range(nlines)]
                    destpath.compare(sorted(expected),path="{}.1_data".format(radix),files_only=True,recursive=False)
                else:
                    radixout = radix
                    
                # Check pymca output (files)
                if cfgfileuse:
                    if addspectra:
                        expected = ["{}_xia{}_{:04d}_0000_{}.edf".format(radixout,det,mapnum,label)\
                                                                        for det in expectedgroups_data\
                                                                        for mapnum in range(nmaps)\
                                                                        for label in fitlabelsfile]
                        expected.extend(["{}_xia{}_{:04d}_0000.cfg".format(radixout,det,mapnum,label)\
                                                                        for det in expectedgroups_data\
                                                                        for mapnum in range(nmaps)])
                    else:
                        expected = ["{}_xia{}_{:04d}_0000_{}.edf".format(radixout,det,mapnum,label)\
                                                                        for det in expectedgroups_data\
                                                                        for mapnum in range(nmaps)\
                                                                        for label in fitlabelsfile]
                        expected.extend(["{}_xia{}_{:04d}_0000.cfg".format(radixout,det,mapnum,label)\
                                                                        for det in expectedgroups_data\
                                                                        for mapnum in range(nmaps)])
                    destpath.compare(sorted(expected),path="{}.1_fit".format(radix),files_only=True,recursive=False)

                # Check top-level output directory (h5 files)
                expected = []
                if cfgfiles:
                    expected.append("{}.1_fit".format(radix))
                if newspectra:
                    expected.append("{}.1_data".format(radix))
                h5file = "{}.h5".format(radix)
                expected.append(h5file)
                destpath.compare(sorted(expected),files_only=True,recursive=False)

                # Check NXprocess groups
                entry = nxprocess.nxentry()
                for name in expected_nxprocess:
                    self.assertTrue(name in entry)
                self.assertEqual(nxprocess.name,expected_nxprocess[-1])

                # Check NXdata groups
                groups, axes = nxresult.regulargriddata(nxprocess)
                self.assertEqual(set(groups.keys()),expectedgroups_result)
                for group,signals in groups.items():
                    if group == 'counters':
                        if quant:
                            expectedsubgroups = counterlabels|calclabels
                        else:
                            expectedsubgroups = counterlabels
                    else:
                        if cfgfileuse:
                            expectedsubgroups = detcounterlabels|fitlabels
                        else:
                            expectedsubgroups = detcounterlabels
                    self.assertEqual(set([sig.name for sig in signals]),expectedsubgroups)

                # Check generated spectra (data)
                if newspectra:
                    # Apply DT correction
                    if dtcorbefore:
                        data0 = data * (stats[...,xiaedf.xiadata.STICR,:]/stats[...,xiaedf.xiadata.STOCR,:])[...,np.newaxis,:]
                    else:
                        data0 = data.copy()
                    
                    # Apply flux normalization
                    if fluxnormbefore:
                        data0 /= ctrs["arr_norm"][...,np.newaxis,np.newaxis]

                    # Add spectra
                    if addspectra:
                        if adddetectorgroups:
                            data0 = np.stack([data0[...,instance.asarray(ind)].sum(axis=-1) for ind in include_detectors],axis=-1)
                        else:
                            data0 = data0[...,alldetectors].sum(axis=-1)[...,np.newaxis]
                    else:
                        data0 = data0[...,alldetectors]
                         
                    # Saved spectra
                    stack = xiaedf.xiastack_radix(os.path.join(destpath.path,"{}.1_data".format(radix)),radixout)
                    data2 = stack.data
                    
                    # Check spectra are equal
                    np.testing.assert_allclose(data0,data2,rtol=1e-6)
                
                # Check fit results
                if cfgfileuse and dtcor:
                    for group,signals in groups.items():
                        if not group.isdetector:
                            continue
                        
                        if group.issum:
                            if adddetectorgroups:
                                if group.number>len(include_detectors):
                                    dets = alldetectors
                                else:
                                    dets = tuple(instance.asarray(include_detectors[group.number-1]).tolist())
                            else:
                                dets = alldetectors 
                        else:
                            dets = (group.number,)
                        logger.debug('Check fit result for sum of detectors {}'.format(dets))
                        info = self.procinfo['detectors'][dets]
                        
                        for signal in signals:
                            if 'xmap' in signal.name:
                                continue
                            dataset = signal.read()
                            if stackdim==1:
                                grpdata = np.moveaxis(dataset,1,0)
                            elif stackdim==2:
                                grpdata = np.moveaxis(dataset,2,0)
                            else:
                                grpdata = dataset[:]

                            self._assert_fitresult(signal.name,grpdata,info)

            #repeat
        #destpath context

    def fitlabels(self,quant=False):
        labels = list(self.procinfo['labels'])
        if quant:
            labels += ['w'+label for label in labels if 'Scatter' not in label]
        return labels

    def qxrfgeometry(self):
        energy = self.procinfo['energy']

        monitor = qxrf.factory("QXRFGeometry",instrument="id21",diodeI0="iodet1",diodeIt="idet",
                                optics="KB",xrfgeometry=None,simplecalibration=True)
        monitor.setreferenceflux(self.procinfo['flux'])
        monitor.setdefaulttime(self.procinfo['time'])
        monitor.diodeI0.gain = 1e8
        monitor.diodeIt.gain = 1e7
        
        info = {"I0_counts":300,"It_counts":30,"time":1,"dark":True,
                "gaindiodeI0":1e8,"gaindiodeIt":1e7}
        monitor.calibrate(**info)
        
        info = {"I0_counts":400000,"It_counts":100000,"time":1,"dark":False,
                "gaindiodeI0":1e8,"gaindiodeIt":1e7,"energy":min(energy)}
        monitor.calibrate(**info)
        
        info = {"I0_counts":200000,"It_counts":100000,"time":1,"dark":False,
                "gaindiodeI0":1e8,"gaindiodeIt":1e7,"energy":max(energy)}
        monitor.calibrate(**info)
        
        return monitor
    
    @contextlib.contextmanager
    def _destpath_context(self):
        destpath = TempDirectory()
        yield destpath
        destpath.cleanup()
      
    def _assert_fitresult(self,grpname,grpdata,info):
        if 'Scatter' in grpname or 'chisq' in grpname:
            return

        grpname = str(grpname)

        m = re.match("Scatter-(Compton|Peak)([0-9]+)",grpname)
        if m:
            grpname = m.group(1)
            if grpname == "Peak":
                grpname = "Rayleigh"
            values1 = [peakareas[0][grpname][int(m.group(2))] for peakareas in info['peakareas']]
            values2 = [peakareas[1][grpname][int(m.group(2))] for peakareas in info['peakareas']]
        else:
            if grpname.startswith("w"):
                grpname = grpname[1:]
                values1 = [massfractions[0][grpname] for massfractions in info['massfractions']]
                values2 = [massfractions[1][grpname] for massfractions in info['massfractions']]
            else:
                values1 = [peakareas[0][grpname] for peakareas in info['peakareas']]
                values2 = [peakareas[1][grpname] for peakareas in info['peakareas']]

        for data,v1,v2 in zip(grpdata,values1,values2):
            mask = data==np.nanmax(data)
            np.testing.assert_allclose(data[~mask],v1,rtol=1e-4)
            np.testing.assert_allclose(data[mask],v2,rtol=1e-4)
            
def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_fluoxas("test_process"))
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
