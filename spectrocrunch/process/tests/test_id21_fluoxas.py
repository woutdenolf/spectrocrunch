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

from ..id21_fluoxas import process
from ...resources import resource_filename
from ...materials import pymca
from ...io import xiaedf

from ...align import types
from ...align.tests import teststack

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

    def gendata(self):
        # Pymca config file
        cfgfile = resource_filename("test/mca.cfg")
        
        # Generate spectra
        spec = np.load(resource_filename("test/mca.npy"))
        spec = np.stack([spec,spec*2,spec*3],axis=1)
        nchan,ndet = spec.shape
        
        data,C,stackdim = teststack.teststack(types.transformationType.translation,nstacks=1)
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

        # Generate data
        path = self.dir.path
        radix = "test"
        stack = xiaedf.xiastack_radix(path,radix)
        xialabels = ["xia{:02d}".format(i) for i in range(ndet)]
        stack.save(data,xialabels,ctrs=ctrs,ctrnames=ctrnames,ctrheaders=ctrheaders)

        return path,radix,cfgfile,data,ctrnames
    
    def test_process(self):
        # Generate data
        sourcepath,scanname,cfgfile,data,ctrnames = self.gendata()
        
        # Raw data 
        nmaps,nlines,nspec,nchan,ndet = data.shape
        scannumbers = [range(nmaps)]
        
        # Process with different settings
        proci = 0
        for alignmethod in [None]:

            alignreference = None
            labels = self.pymcaoutlabels(cfgfile)
            for label in labels:
                if "_K" in label:
                    alignreference = label.replace("_","-")
            refimageindex = 0
            exclude_detectors = [1]
            
            outpath = "crunched{}".format(proci)
            destpath = os.path.join(sourcepath,outpath)
            
            for skippre in [False,True]:
                process(sourcepath,destpath,scanname,scannumbers,cfgfile,\
                        alignmethod=alignmethod,alignreference=alignreference,\
                        refimageindex=refimageindex,dtcor=False,plot=False,\
                        exclude_detectors=exclude_detectors,skippre=skippre)
                
                # Check generated spectra (files)
                expected = ["{}_raw_xiaS1_{:04d}_0000_{:04d}.edf".format(scanname,mapnum,linenum) for mapnum in range(nmaps) for linenum in range(nlines)]
                self.dir.compare(sorted(expected),path=(outpath,"{}_data".format(scanname)),files_only=True,recursive=False)

                # Check pymca output (files)
                expected = ["{}_raw_xiaS1_{:04d}_0000_{}.edf".format(scanname,mapnum,label) for mapnum in range(nmaps) for label in labels]
                self.dir.compare(sorted(expected),path=(outpath,"{}_fit".format(scanname)),files_only=True,recursive=False)
                
                # Check hdf5 output (files)
                expected = ["{}_data".format(scanname),"{}_fit".format(scanname)]
                expected += ["{}.h5".format(scanname),"{}.json".format(scanname)]
                if alignmethod is not None and alignreference is not None:
                    expected.append("{}.align.h5".format(scanname))
                    expected.append("{}.replace.h5".format(scanname))
                self.dir.compare(sorted(expected),path=outpath,files_only=True,recursive=False)

                # Check generated spectra (data)
                stack = xiaedf.xiastack_radix(os.path.join(destpath,"{}_data".format(scanname)),"{}_raw".format(scanname))
                data2 = stack.data[...,0]
                np.testing.assert_array_equal(data[...,0]+data[...,2],data2)
                
            proci += 1
        

def test_suite_all():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    #testSuite.addTest(test_fluoxas("test_process"))
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite_all()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
