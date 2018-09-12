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
from testfixtures import TempDirectory
from PyMca5.PyMcaIO import ConfigDict
import numpy as np
import os
import contextlib
import itertools
import logging
import re

from ..fluoxas import process
from ...io import xiaedf
from ...io import nexus
from ...align import types
from ...common import instance
from ...geometries import qxrf
from ...h5stacks import get_hdf5_imagestacks
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

logger = logging.getLogger(__loader__.fullname)

class test_fluoxas(unittest.TestCase):

    def setUp(self):
        self.dir = TempDirectory()
    
    def create_procinfo(self):
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
        pymcahandle = pymca.PymcaHandle(energy=8.0,flux=1e9,time=0.1,snip=False,continuum=1,escape=False,linear=True,noisepropagation=False)

        ndet = 3
        detectorsinsum = [0,2]
        samult = [1,0.8,0.5]
        samult.append(sum(i for i in samult))
        energy = [8,8.5,9]
        
        detectornames = ['detector{}'.format(i) for i in range(ndet)+['sum']]
        detectors = []
        for det,sam in zip(detectornames,samult):
            geometry = xrfgeometries.factory("sxm120",detector=leia,
                                         source=synchrotron,detectorposition=-15.)
            geometry.solidangle *= sam

            cfgfile = os.path.join(self.dir.path,det+'.cfg')
            detector = {'cfgfile':cfgfile,'geometry':geometry,'mca':[],'peakareas':[],'massfractions':[]}
            detectors.append(detector)

            for en in energy:
                pymcahandle.energy = en
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
                    
                    pymcahandle.loadfrompymca(cfgfile)
                    pymcahandle.savepymca(cfgfile+'_')
                    from ...materials.pymcadiff import diff
                    diff(cfgfile,cfgfile+'_')
                    #exit()
                    
                    mca = pymcahandle.mca()+1
                    mcas.append(mca)

                    # This check is done more rigorously in test_pymca
                    fpeakareas = pymcahandle.xraygroupareas()
                    fmassfractions = pymcahandle.sample.elemental_massfractions()
                    pymcahandle.setdata(mca)
                    fitresult = pymcahandle.fit(loadfromfit=False)
                    for k,v in fitresult["fitareas"].items():
                        np.testing.assert_allclose(v,fpeakareas[k],rtol=1e-2)
                    for k,v in fitresult["massfractions"].items():
                        np.testing.assert_allclose(v,fmassfractions[k.element],rtol=2e-2)
                    peakareas.append(fitresult["fitareas"])
                    massfractions.append(fitresult["massfractions"])
        exit()
        self.procinfo = {}
        self.procinfo['include_detectors'] = ([2],detectorsinsum)
        self.procinfo['flux'] = pymcahandle.flux
        self.procinfo['time'] = pymcahandle.time
        self.procinfo['detectorsum'] = detectors.pop(-1)
        self.procinfo['detectors'] = detectors
        self.procinfo['energy'] = energy
        
    def tearDown(self):
        pass#self.dir.cleanup()

    def fitlabelsfile(self,quant=False):
        config = ConfigDict.ConfigDict()
        config.read(self.procinfo['detectorsum']['cfgfile'])

        labels = ["{}_{}".format(e,line) for e,lines in config["peaks"].items() for line in lines]
        if quant:
            labels = labels + ["w{}_{}".format(e,line) for e,lines in config["peaks"].items() for line in lines]
        
        if config["fit"]["scatterflag"]:
            n = sum(instance.asarray(config["fit"]["energyflag"]) &\
                   instance.asarray(config["fit"]["energyscatter"]))
            labels += ["Scatter_Compton{:03d}".format(i) for i in range(n)]
            labels += ["Scatter_Peak{:03d}".format(i) for i in range(n)]

        return labels

    def qxrfgeometry(self):
        energy = self.procinfo['energy']

        monitor = qxrf.factory("QXRFGeometry",instrument="id21",diodeI0="iodet1",diodeIt="idet",
                                optics="KB",xrfgeometry=None,simplecalibration=True)
        monitor.setreferenceflux(self.procinfo['flux'])
        monitor.setdefaulttime(self.procinfo['time'])
        monitor.diodeI0.gain = 1e8
        monitor.diodeIt.gain = 1e7
        
        info = {"I0_counts":300,"It_counts":30,"time":1,"dark":True,"gaindiodeI0":1e8,"gaindiodeIt":1e7}
        monitor.calibrate(**info)
        
        info = {"I0_counts":400000,"It_counts":100000,"time":1,"dark":False,"gaindiodeI0":1e8,"gaindiodeIt":1e7,"energy":min(energy)}
        monitor.calibrate(**info)
        
        info = {"I0_counts":200000,"It_counts":100000,"time":1,"dark":False,"gaindiodeI0":1e8,"gaindiodeIt":1e7,"energy":max(energy)}
        monitor.calibrate(**info)
        
        return monitor
    
    def gendata(self):
        self.create_procinfo()
        qxrfgeometry = self.qxrfgeometry()
        refflux = qxrfgeometry.reference.to("hertz").magnitude
        expotime = qxrfgeometry.defaultexpotime.to("seconds").magnitude

        # 3 maps of a moving hotspot
        nlines,nspec = 7,6
        nmaps = len(self.procinfo['energy'])
        ndet = len(self.procinfo['detectors'])
        nchan = len(self.procinfo['detectorsum']['mca'][0][0])
        data = np.ones((nmaps,nlines,nspec,nchan,ndet))
        off = max(min(nlines,nspec)-nmaps,0)//2
        for idet,det in enumerate(self.procinfo['detectors']):
            for imap,spectra in enumerate(det['mca']):
                spec1,spec2 = spectra
                data[imap,...,idet] = spec1
                data[imap,imap+off,imap+off,:,idet] = spec2

        # Generate counter headers
        energy = self.procinfo['energy']
        ctrheaders = np.vectorize(lambda e:{"DCM_Energy":e,"time":expotime},otypes=[object])(energy)
        
        # Init counters
        ctrs = {}
        
        # Apply flux decay
        rflux = np.linspace(1,0.5,nmaps*nlines*nspec).reshape((nmaps,nlines,nspec))
        flux = rflux*refflux
        data *= rflux[...,np.newaxis,np.newaxis]

        ctrs["arr_iodet"] = np.ones((nmaps,nlines,nspec))
        ctrs["arr_idet"] = np.ones((nmaps,nlines,nspec))
        ctrs["arr_norm"] = np.ones((nmaps,nlines,nspec))
        for i,en in enumerate(energy):
            ctrs["arr_iodet"][i,...] = qxrfgeometry.diodeI0.fluxtocps(en,flux[i,...])*expotime
            ctrs["arr_idet"][i,...] = qxrfgeometry.diodeIt.fluxtocps(en,flux[i,...])*expotime
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
            OCR = ICR*(1-0.2*i/(ndet-1.))# DT = 10*i %
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

        return path,radix,data,stats,ctrs,qxrfgeometry
    
    @contextlib.contextmanager
    def env_destpath(self):
        self.destpath = TempDirectory()
        yield
        self.destpath.cleanup()
      
    def _assert_fitresult(self,grpname,grpdata,info):
        if 'Scatter' in grpname:
            return

        grpname = str(grpname)
        print grpname
        
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
            print data
            print v1,v2
            mask = data==np.nanmax(data)
            np.testing.assert_allclose(data[~mask],v1,rtol=1e-3)
            np.testing.assert_allclose(data[mask],v2,rtol=1e-3)
                                      
    def test_process(self):
        # Generate data
        sourcepath,radix,data,stats,ctrs,qxrfgeometry = self.gendata()
        
        # Raw data 
        nmaps,nlines,nspec,nchan,ndet = data.shape
        scannumbers = [range(nmaps)]

        parameters = [(None,"max"),(True,False),self.procinfo['include_detectors'],(True,False),
                      (True,False),(True,False),(True,False),(2,),(True,False)]
        for combination in itertools.product(*parameters):
            alignmethod,cfgfileuse,include_detectors,adddetectors,\
            addbeforefit,quant,dtcor,stackdim,correctspectra = combination
            if not cfgfileuse and alignmethod is not None:
                continue
            adddetectors = len(include_detectors)>1 and adddetectors
            addbefore = adddetectors and addbeforefit
            fluxnormbefore = quant and correctspectra
            dtcorbefore = dtcor and (correctspectra or addbefore)
            newspectra = addbefore or fluxnormbefore or dtcorbefore
                
            if quant:
                geom = qxrfgeometry
                if addbefore:
                    geom.xrfgeometries = self.procinfo['detectorsum']['geometry']
                else:
                    geom.xrfgeometries = [self.procinfo['detectors'][i]['geometry'] for i in include_detectors]
                prealignnormcounter = None
            else:
                geom = None
                prealignnormcounter = "arr_norm"

            if cfgfileuse:
                if addbefore:
                    cfgfiles = self.procinfo['detectorsum']['cfgfile']
                else:
                    cfgfiles = [self.procinfo['detectors'][i]['cfgfile'] for i in include_detectors]
            else:
                cfgfiles = None

            alignreference = None
            fitlabelsfile = self.fitlabelsfile(quant=quant)
            fitlabels = set([label.replace("_","-") for label in fitlabelsfile])
            detcounterlabels = set(['xmap_icr','xmap_ocr','xmap_x1c','xmap_x2c'])
            counterlabels = set(['arr_iodet','arr_idet','arr_norm'])
            calclabels = set(['calc_transmission','calc_absorbance','calc_flux0','calc_fluxt'])
            for label in fitlabels:
                if not 'Scatter' in label:
                    alignreference = label
                    break
            refimageindex = 0

            if adddetectors:
                alignref = "/detectorsum/{}".format(alignreference)
            else:
                alignref = "/detector{}/{}".format(include_detectors[0],alignreference)
            
            with self.env_destpath():
                for skippre in [False,]:
                    logger.debug("skippre = {}".format(skippre))
                    logger.debug("stackdim = {}".format(stackdim))
                    logger.debug("cfgfiles = {}".format(cfgfiles))
                    logger.debug("dtcor = {}".format(dtcor))
                    logger.debug("adddetectors = {}".format(adddetectors))
                    logger.debug("addbeforefit = {}".format(addbeforefit))
                    logger.debug("include_detectors = {}".format(include_detectors))
                    logger.debug("quant = {}".format(quant))
                    logger.debug("correctspectra = {}".format(correctspectra))
                    logger.debug("newspectra = {}".format(newspectra))
                    
                    parameters = {}
                    parameters["alignmethod"] = alignmethod
                    parameters["alignreference"] = alignref
                    parameters["refimageindex"] = refimageindex
                    parameters["dtcor"] = dtcor
                    parameters["plot"] = False
                    parameters["adddetectors"] = adddetectors
                    parameters["addbeforefit"] = addbeforefit
                    parameters["qxrfgeometry"] = geom
                    parameters["correctspectra"] = correctspectra
                    parameters["replacenan"] = bool(alignmethod)
                    parameters["prealignnormcounter"] = prealignnormcounter
                    parameters["stackdim"] = stackdim
                    parameters["include_detectors"] = include_detectors
                    parameters["skippre"] = skippre
                    parameters["instrument"] = "id21"
                    parameters["counters"] = ["arr_norm"]

                    process(sourcepath,self.destpath.path,radix,scannumbers,cfgfiles,**parameters)

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

                        if addbefore:
                            expected = ["{}_xiaS1_{:04d}_0000_{:04d}.edf".format(radixout,mapnum,linenum)\
                                                                        for mapnum in range(nmaps)\
                                                                        for linenum in range(nlines)]
                        else:
                            expected = ["{}_xia{:02d}_{:04d}_0000_{:04d}.edf".format(radixout,det,mapnum,linenum)\
                                                                        for det in include_detectors\
                                                                        for mapnum in range(nmaps)\
                                                                        for linenum in range(nlines)]
                        self.destpath.compare(sorted(expected),path="{}_data".format(radix),files_only=True,recursive=False)
                    else:
                        radixout = radix
                        
                    # Check pymca output (files)
                    if cfgfileuse:
                        if addbefore:
                            expected = ["{}_xiaS1_{:04d}_0000_{}.edf".format(radixout,mapnum,label)\
                                                                            for mapnum in range(nmaps)\
                                                                            for label in fitlabelsfile]
                            expected.extend(["{}_xiaS1_{:04d}_0000.cfg".format(radixout,mapnum,label)\
                                                                            for mapnum in range(nmaps)])
                        else:
                            expected = ["{}_xia{:02d}_{:04d}_0000_{}.edf".format(radixout,det,mapnum,label)\
                                                                            for det in include_detectors\
                                                                            for mapnum in range(nmaps)\
                                                                            for label in fitlabelsfile]
                            expected.extend(["{}_xia{:02d}_{:04d}_0000.cfg".format(radixout,det,mapnum,label)\
                                                                            for det in include_detectors\
                                                                            for mapnum in range(nmaps)])
                        self.destpath.compare(sorted(expected),path="{}_fit".format(radix),files_only=True,recursive=False)

                    # Check top-level output directory (h5 files)
                    expected = []
                    if cfgfiles:
                        expected.append("{}_fit".format(radix))
                    if newspectra:
                        expected.append("{}_data".format(radix))
                    h5file = "{}.h5".format(radix)
                    expected.append(h5file)
                    expected.append("{}.json".format(radix))
                    if prealignnormcounter is not None and cfgfileuse:
                        h5file = "{}.norm.h5".format(radix)
                        expected.append(h5file)
                    if alignmethod is not None and alignreference is not None:
                        expected.append("{}.align.h5".format(radix))
                        h5file = "{}.replace.h5".format(radix)
                        expected.append(h5file)
                    self.destpath.compare(sorted(expected),files_only=True,recursive=False)

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
                        data0 = data0[...,include_detectors]
                        if addbefore:
                            data0 = data0.sum(axis=-1)[...,np.newaxis]

                        # Saved spectra
                        stack = xiaedf.xiastack_radix(os.path.join(self.destpath.path,"{}_data".format(radix)),radixout)
                        data2 = stack.data
                        
                        # Check spectra are equal
                        np.testing.assert_allclose(data0,data2,rtol=1e-6)
                    
                    # Check HDF5 groups
                    h5file = os.path.join(self.destpath.path,h5file)
                    stacks, axes, procinfo = get_hdf5_imagestacks.get_hdf5_imagestacks(h5file,['^(counters|(detector([0-9]+|sum)))$'])
                    
                    expectedgroups = ['counters']
                    if adddetectors:
                        expectedgroups.append('detectorsum')
                    else:
                        expectedgroups += ['detector{}'.format(i) for i in include_detectors]
                    expectedgroups = set(expectedgroups)
                    self.assertEqual(set(stacks.keys()),expectedgroups)

                    for detector,stack in stacks.items():
                        if detector == 'counters':
                            if quant:
                                expectedgroups = counterlabels|calclabels
                            else:
                                expectedgroups = counterlabels
                        else:
                            if cfgfileuse:
                                expectedgroups = detcounterlabels|fitlabels
                            else:
                                expectedgroups = detcounterlabels
                        self.assertEqual(set(stack.keys()),expectedgroups)
                        
                    # Check fit results
                    if cfgfileuse and dtcor:
                        with nexus.File(h5file,mode='r') as f:
                            for detector,stack in stacks.items():
                                if detector == 'detectorsum':
                                    info = self.procinfo['detectorsum']
                                else:
                                    info = self.procinfo['detectors'][int(detector[8:])]

                                for grpname,grp in stack.items():
                                    if 'xmap' in grp:
                                        continue
                                    dataset,_,_ = nexus.parse_NXdata(f[grp])
                                    if stackdim==1:
                                        grpdata = np.moveaxis(dataset,1,0)
                                    elif stackdim==2:
                                        grpdata = np.moveaxis(dataset,2,0)
                                    else:
                                        grpdata = dataset[:]
       
                                    self._assert_fitresult(grpname,grpdata,info)
                                
                    continue
                    if cfgfileuse:
                        h5file = os.path.join(self.destpath.path,h5file)
                        stacks, axes, procinfo = get_hdf5_imagestacks.get_hdf5_imagestacks(h5file,["detectorsum"])
                        data4 = None
                        with nexus.File(h5file,mode='r') as f:
                            for stack in stacks.values():
                                for grp in stack.values():
                                    #if "xmap" in grp:
                                    #    continue
                                    dataset,_,_ = nexus.parse_NXdata(f[grp])
                                    if stackdim==1:
                                        data3 = np.moveaxis(dataset,1,0)
                                    elif stackdim==2:
                                        data3 = np.moveaxis(dataset,2,0)
                                    else:
                                        data3 = dataset[:]
                                    if data4 is None:
                                        data4 = data3
                                    else:
                                        r = data4/data3 # ratio's of elements
                                        
                                        #TODO: large differences in ratio's, why?

                                        try:
                                            np.testing.assert_allclose(np.nanmin(r),np.nanmax(r),rtol=1e-6)
                                        except Exception as e:
                                            print h5file
                                            print grp
                                            print (np.nanmin(r)-np.nanmax(r))/np.nanmax(r)
                                            for i in range(nmaps):
                                                print data4[i,...]
                                                print data3[i,...]
                                                print r[i,...]
                                            raise e
                    # Finished one condition set
                                
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
