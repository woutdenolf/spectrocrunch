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

from ..id21_fluoxas import process
from ...io import xiaedf
from ...io import nexus
from ...align import types
from ...common import instance
from ..id21_quant import FluxMonitor
from ...h5stacks.get_hdf5_imagestacks import get_hdf5_imagestacks as getstacks
from ...materials.tests.xrf_setup import pymcahandle

class test_fluoxas(unittest.TestCase):

    def setUp(self):
        self.dir = TempDirectory()
        self.cfgfile = os.path.join(self.dir.path,"mca.cfg")
        pymcahandle.savepymca(self.cfgfile)

    def tearDown(self):
        self.dir.cleanup()

    def pymcaoutlabels(self,quant=False):
        config = ConfigDict.ConfigDict()
        config.read(self.cfgfile)

        labels = ["{}_{}".format(e,line) for e,lines in config["peaks"].items() for line in lines]
        if quant:
            labels = labels + ["C({}_{})".format(e,line) for e,lines in config["peaks"].items() for line in lines]
        
        if config["fit"]["scatterflag"]:
            n = sum(instance.asarray(config["fit"]["energyflag"]) &\
                   instance.asarray(config["fit"]["energyscatter"]))
            labels += ["Scatter_Compton{:03d}".format(i) for i in range(n)]
            labels += ["Scatter_Peak{:03d}".format(i) for i in range(n)]

        return labels

    def pymcagetenergy(self):
        config = ConfigDict.ConfigDict()
        config.read(self.cfgfile)
        return float(instance.asarray(config["fit"]["energy"])[0])

    def fluxmonitor(self):
        energy = self.pymcagetenergy()
        monitor = FluxMonitor(iodetname="iodet1",focussed=True,xrfdetector="leia",xrfgeometry="sxm120")
        monitor.setdark(300,None,gainiodet=1e8)
        monitor.setcalib(energy-5,0.5,gainiodet=1e8)
        monitor.setcalib(energy+5,0.5,gainiodet=1e8)
        monitor.setreferenceflux(1e9)
        monitor.settime(0.1)
        return monitor
    
    def gendata(self):
        # Generate spectra
        spec =  pymcahandle.mca()+1
        spec = np.stack([spec,spec*2,spec*3],axis=1)
        nchan,ndet = spec.shape
        
        nmaps,nlines,nspec = 3,7,6
        data = np.ones((nmaps,nlines,nspec))
        for i in range(nmaps):
            t = i-nmaps//2
            data[i,t,t] = 2
        data = np.outer(data,spec).reshape(nmaps,nlines,nspec,nchan,ndet)

        # Generate counter headers
        energy = self.pymcagetenergy()
        energy = np.linspace(energy,energy+0.001,nmaps)
        ctrheaders = np.vectorize(lambda e:{"DCM_Energy":e},otypes=[object])(energy)
        
        # Init counters
        ctrs = {}

        # Apply flux
        fluxmonitor = self.fluxmonitor()
        refflux = fluxmonitor.reference.to("hertz").magnitude
        flux = np.linspace(refflux,refflux*0.5,nmaps*nlines*nspec).reshape((nmaps,nlines,nspec))
        data *= flux[...,np.newaxis,np.newaxis]
        
        ctrs["arr_iodet"] = np.ones((nmaps,nlines,nspec))
        ctrs["arr_idet"] = np.ones((nmaps,nlines,nspec))
        ctrs["arr_fdet"] = np.ones((nmaps,nlines,nspec))
        for i,en in enumerate(energy):
            ctrs["arr_iodet"][i,...] = fluxmonitor.fluxtocps(en,flux[i,...])*fluxmonitor.defaulttime.to("seconds").magnitude
            ctrs["arr_fdet"][i,...] = flux[i,...]/refflux
            op,fref,tref,traw = fluxmonitor.xrfnormop(en)
            self.assertEqual(fref,refflux)
            np.testing.assert_allclose(ctrs["arr_fdet"][i,...],op(ctrs["arr_iodet"][i,...]))

        # Deadtime corrected counters
        a = nchan/2
        b = a+100
        for i in range(ndet):
            ctrs["xmap_x1c_{:02d}".format(i)] = data[...,a:b,i].sum(axis=-1)

        b = a
        a = a-100
        for i in range(ndet):
            ctrs["xmap_x2c_{:02d}".format(i)] = data[...,a:b,i].sum(axis=-1)

        # Transmission
        ctrs["arr_idet"] = np.zeros((nmaps,nlines,nspec)) # no transmission

        # Apply deadtime
        stats = np.zeros((nmaps,nlines,nspec,xiaedf.xiadata.NSTATS,ndet),dtype=data.dtype)
        for j in range(ndet):
            ICR = data[...,j].sum(axis=-1)
            OCR = ICR*(1-j/50.)# DT = 2*j %
            
            ctrs["xmap_icr_{:02d}".format(i)] = ICR
            ctrs["xmap_ocr_{:02d}".format(i)] = OCR
            
            data[...,j] = data[...,j]*(OCR/ICR)[...,np.newaxis]

            stats[...,xiaedf.xiadata.STDET,j] = j
            stats[...,xiaedf.xiadata.STEVT,j] = ICR # % Not sure
            stats[...,xiaedf.xiadata.STICR,j] = ICR
            stats[...,xiaedf.xiadata.STOCR,j] = OCR
            stats[...,xiaedf.xiadata.STDT,j] = 100-OCR*100./ICR # %
            stats[...,xiaedf.xiadata.STLT,j] = OCR*1000./ICR # 1000 msec RT

        # Generate data
        path = self.dir.path
        radix = "test"
        stack = xiaedf.xiastack_radix(path,radix)
        xialabels = ["xia{:02d}".format(i) for i in range(ndet)]
        stack.save(data,xialabels,stats=stats,ctrs=np.stack(ctrs.values(),axis=-1),ctrnames=ctrs.keys(),ctrheaders=ctrheaders)

        return path,radix,data,stats,ctrs,fluxmonitor
    
    @contextlib.contextmanager
    def env_destpath(self):
        self.destpath = TempDirectory()
        yield
        self.destpath.cleanup()

    def test_process(self):
        
        # Generate data
        sourcepath,radix,data,stats,ctrs,fluxmonitor = self.gendata()
        
        # Raw data 
        nmaps,nlines,nspec,nchan,ndet = data.shape
        scannumbers = [range(nmaps)]
        
        # Fixed process parameters
        alignreference = None
        labels = self.pymcaoutlabels()
        for label in labels:
            if "_K" in label:
                alignreference = label.replace("_","-")
        refimageindex = 0

        # Process with different settings
        for alignmethod in ["max",None]:
        
            for cfgfileuse in [None,self.cfgfile]:
                if cfgfileuse is None and alignmethod is not None:
                    continue

                for include_detectors in [[2],[0,2]]:
                
                    for addbeforefit in [True,False]:
                    
                        for quant in [True,False]:
                            if quant:
                                monitor = fluxmonitor
                                prealignnormcounter = None
                            else:
                                monitor = None
                                prealignnormcounter = "arr_fdet"
                            
                            for dtcor in [True,False]:
                                dtcor_onspectra = dtcor and len(include_detectors)>1

                                if dtcor_onspectra:
                                    radixout = "{}_{}".format(radix,"dtcor")
                                else:
                                    radixout = radix

                                for stackdim in [2,1,0]:
                                
                                    with self.env_destpath():
                                        for skippre in [False]:
                                            addbeforefit_onspectra = addbeforefit and len(include_detectors)>1
                                            
                                            #print "alignmethod=",alignmethod
                                            #print "addbeforefit=",addbeforefit
                                            #print "addbeforefit_onspectra=",addbeforefit_onspectra
                                            #print "dtcor=",dtcor
                                            #print "dtcor_onspectra=",dtcor_onspectra
                                            #print "skippre=",skippre
                                            #print "quant=",quant
                                            #print "stackdim=",stackdim
                                            
                                            process(sourcepath,self.destpath.path,radix,scannumbers,cfgfileuse,\
                                                    alignmethod=alignmethod,alignreference=alignreference,\
                                                    refimageindex=refimageindex,dtcor=dtcor,plot=False,\
                                                    addbeforefit=addbeforefit,fluxmonitor=monitor,replacenan=True,\
                                                    prealignnormcounter=prealignnormcounter,stackdim=stackdim,\
                                                    include_detectors=include_detectors,skippre=skippre)

                                            # Check generated spectra (files)
                                            newspectra = dtcor_onspectra or addbeforefit_onspectra or quant
                                            
                                            if newspectra:
                                                if addbeforefit_onspectra:
                                                    expected = ["{}_xiaS1_{:04d}_0000_{:04d}.edf".format(radixout,mapnum,linenum)\
                                                                                                for mapnum in range(nmaps)\
                                                                                                for linenum in range(nlines)]
                                                else:
                                                    expected = ["{}_xia{:02d}_{:04d}_0000_{:04d}.edf".format(radixout,det,mapnum,linenum)\
                                                                                                for det in include_detectors\
                                                                                                for mapnum in range(nmaps)\
                                                                                                for linenum in range(nlines)]
                                                self.destpath.compare(sorted(expected),path="{}_data".format(radix),files_only=True,recursive=False)
                                                
                                            # Check pymca output (files)
                                            if cfgfileuse is not None:
                                                labels = self.pymcaoutlabels(quant=quant)
                                                if addbeforefit_onspectra:
                                                    expected = ["{}_xiaS1_{:04d}_0000_{}.edf".format(radixout,mapnum,label)\
                                                                                                    for mapnum in range(nmaps)\
                                                                                                    for label in labels]
                                                else:
                                                    expected = ["{}_xia{:02d}_{:04d}_0000_{}.edf".format(radixout,det,mapnum,label)\
                                                                                                    for det in include_detectors\
                                                                                                    for mapnum in range(nmaps)\
                                                                                                    for label in labels]
                                                self.destpath.compare(sorted(expected),path="{}_fit".format(radix),files_only=True,recursive=False)
                                            
                                            # Check top-level output directory (files)
                                            expected = []
                                            if cfgfileuse is not None:
                                                expected.append("{}_fit".format(radix))
                                            if newspectra:
                                                expected.append("{}_data".format(radix))
                                            h5file = "{}.h5".format(radix)
                                            expected.append(h5file)
                                            expected.append("{}.json".format(radix))
                                            if prealignnormcounter is not None:
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
                                                if dtcor_onspectra:
                                                    data0 = data * (stats[...,xiaedf.xiadata.STICR,:]/stats[...,xiaedf.xiadata.STOCR,:])[...,np.newaxis,:]
                                                else:
                                                    data0 = data.copy()
                                                
                                                # Apply flux normalization
                                                if quant:
                                                    data0 /= ctrs["arr_fdet"][...,np.newaxis,np.newaxis]

                                                # Add spectra
                                                if addbeforefit_onspectra:
                                                    data0 = data0[...,0]+data0[...,2]
                                                    data0 = data0[...,np.newaxis]
                                                else:
                                                    data0 = data0[...,include_detectors]

                                                # Saved spectra
                                                stack = xiaedf.xiastack_radix(os.path.join(self.destpath.path,"{}_data".format(radix)),radixout)
                                                data2 = stack.data
                                                
                                                # Check spectra are equal
                                                np.testing.assert_allclose(data0,data2,rtol=1e-6)

                                            # Check element ratio's the same in all pixels
                                            if cfgfileuse is not None:
                                                h5file = os.path.join(self.destpath.path,h5file)
                                                stacks, axes = getstacks(h5file,["detectorsum"])
                                                data4 = None
                                                with nexus.File(h5file,mode='r') as f:
                                                    for stack in stacks.values():
                                                        for grp in stack.values():
                                                            if "xmap" in grp:
                                                                continue
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
                                                                
                                                                #for i in range(nmaps):
                                                                #    v = np.median(r[i,...])
                                                                #    ind = np.unravel_index(np.argmax(np.abs(r[i,...]-v)), r.shape[1:])
                                                                #    r[i,ind[0],ind[1]]= np.nan
                                                                
                                                                #try:
                                                                np.testing.assert_allclose(np.nanmin(r),np.nanmax(r),rtol=1e-1)
                                                                #except Exception as e:
                                                                #    print h5file
                                                                #    print grp
                                                                #    for i in range(nmaps):
                                                                #        print r[i,...]
                                                                #    raise e
       
                                
def test_suite_all():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_fluoxas("test_process"))
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite_all()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
