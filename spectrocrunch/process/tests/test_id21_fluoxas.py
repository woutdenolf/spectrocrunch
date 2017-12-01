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
from ...resources import resource_filename
from ...materials import pymca
from ...io import xiaedf
from ...align import types
from ...align.tests import teststack
from ..id21_quant import FluxMonitor

class test_fluoxas(unittest.TestCase):

    def setUp(self):
        self.dir = TempDirectory()

    def tearDown(self):
        self.dir.cleanup()
    
    def pymcaoutlabels(self,cfgfile):
        config = ConfigDict.ConfigDict()
        config.read(cfgfile)

        labels = ["{}_{}".format(e,line) for e,lines in config["peaks"].items() for line in lines]
        
        #if config["fit"]["scatterflag"]:
        #    labels += ["Scatter_Compton{:03d}".format(i) for i,b in enumerate(config["fit"]["energyflag"]) if b]
        #    labels += ["Scatter_Peak{:03d}".format(i) for i,b in enumerate(config["fit"]["energyflag"]) if b]
        labels += ["Scatter_Compton000","Scatter_Peak000"]
        
        return labels

    def pymcagetenergy(self,cfgfile):
        config = ConfigDict.ConfigDict()
        config.read(cfgfile)
        
        return config["fit"]["energy"][0]

    def fluxmonitor(self,cfgfile):
        energy = self.pymcagetenergy(cfgfile)
        monitor = FluxMonitor()
        monitor.setdiode("iodet1")
        monitor.manualdark(300,1e5)
        monitor.manualcalib(energy-5,0.5,1e5)
        monitor.manualcalib(energy+5,0.5,1e5)
        return monitor
        
    def gendata(self):
        # Pymca config file
        cfgfile = resource_filename("test/mca.cfg")
        
        # Generate spectra
        spec = np.load(resource_filename("test/mca.npy"))
        spec = np.stack([spec,spec*2,spec*3],axis=1)
        nchan,ndet = spec.shape
        
        data,C,stackdim = teststack.teststack(types.transformationType.translation,\
                                        nstacks=1,ndim1 = 11,ndim2 = 15,nimages = 3)

        data = data[0]
        if stackdim==1:
            data = np.moveaxis(data,1,0)
        elif stackdim==2:
            data = np.moveaxis(data,2,0)
        nmaps,nlines,nspec = data.shape
        data = np.outer(data,spec).reshape(nmaps,nlines,nspec,nchan,ndet)

        # Generate counters
        ctrnames = ["xmap_x1_00","xmap_x1_01","xmap_x1_02","arr_iodet","arr_idet"]
        ncounters = len(ctrnames)
        ctrs = np.ones((nmaps,nlines,nspec,ncounters))
        a = nchan/2
        b = a+100
        for i in range(ndet):
            ctrs[...,i] = data[...,a:b,i].sum(axis=-1)

        # Generate counter headers
        energy = self.pymcagetenergy(cfgfile)
        energy = np.linspace(energy,energy+0.001,nmaps)
        ctrheaders = np.vectorize(lambda e:{"DCM_Energy":e},otypes=[object])(energy)

        # Generate statistics
        stats = np.zeros((nmaps,nlines,nspec,xiaedf.xiadata.NSTATS,ndet),dtype=data.dtype)
        for j in range(ndet):
            ICR = data[...,j].sum(axis=-1)
            OCR = ICR*(1-j/50.)# DT = 2*j %

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
        stack.save(data,xialabels,stats=stats,ctrs=ctrs,ctrnames=ctrnames,ctrheaders=ctrheaders)

        return path,radix,cfgfile,data,stats,ctrs,ctrnames
    
    @contextlib.contextmanager
    def env_destpath(self):
        self.destpath = TempDirectory()
        yield
        self.destpath.cleanup()

    def test_process(self):
        # Generate data
        sourcepath,radix,cfgfile,data,stats,ctrs,ctrnames = self.gendata()
        
        # Raw data 
        nmaps,nlines,nspec,nchan,ndet = data.shape
        scannumbers = [range(nmaps)]
        
        # Fixed process parameters
        alignreference = None
        labels = self.pymcaoutlabels(cfgfile)
        for label in labels:
            if "_K" in label:
                alignreference = label.replace("_","-")
        refimageindex = 0
        exclude_detectors = [1]
        include_detectors = [i for i in range(ndet) if i not in exclude_detectors]
        
        # Process with different settings
        for alignmethod in [None]:
            
            for addbeforefit in [True,False]:
            
                for quant in [False]:
                    if quant:
                        fluxmonitor = self.fluxmonitor(cfgfile)
                    else:
                        fluxmonitor = None
                        
                    for dtcor in [False,True]:
                        if dtcor:
                            radixout = "{}_{}".format(radix,"dtcor")
                        else:
                            radixout = radix

                        with self.env_destpath():
                            for skippre in [False,True]:
                                print "alignmethod=",alignmethod
                                print "addbeforefit=",addbeforefit
                                print "dtcor=",dtcor
                                print "skippre=",skippre
                                print "quant=",quant
                                
                                process(sourcepath,self.destpath.path,radix,scannumbers,cfgfile,\
                                        alignmethod=alignmethod,alignreference=alignreference,\
                                        refimageindex=refimageindex,dtcor=dtcor,plot=True,\
                                        addbeforefit=addbeforefit,fluxmonitor=fluxmonitor,\
                                        exclude_detectors=exclude_detectors,skippre=skippre)
                                
                                # Check generated spectra (files)
                                newspectra = dtcor or addbeforefit
                                
                                if newspectra:
                                    if addbeforefit:
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
                                if addbeforefit:
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
                                expected = ["{}_fit".format(radix)]
                                if newspectra:
                                    expected += ["{}_data".format(radix)]
                                expected += ["{}.h5".format(radix),"{}.json".format(radix)]
                                if alignmethod is not None and alignreference is not None:
                                    expected.append("{}.align.h5".format(radix))
                                    expected.append("{}.replace.h5".format(radix))
                                self.destpath.compare(sorted(expected),files_only=True,recursive=False)

                                # Check generated spectra (data)
                                if newspectra:
                                    # Saved spectra
                                    stack = xiaedf.xiastack_radix(os.path.join(self.destpath.path,"{}_data".format(radix)),radixout)
                                    data2 = stack.data
                                        
                                    # Apply DT correction
                                    if dtcor:
                                        data0 = data * (stats[...,xiaedf.xiadata.STICR,:]/stats[...,xiaedf.xiadata.STOCR,:])[...,np.newaxis,:]
                                    else:
                                        data0 = data
                                    
                                    # Add spectra
                                    if addbeforefit:
                                        data0 = data0[...,0]+data0[...,2]
                                        data0 = data0[...,np.newaxis]
                                    else:
                                        data0 = data0[...,include_detectors]

                                    # Check spectra are equal
                                    np.testing.assert_array_equal(data0,data2)
                                
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
