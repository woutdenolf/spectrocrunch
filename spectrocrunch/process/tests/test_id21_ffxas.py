# -*- coding: utf-8 -*-
#
#   Copyright (C) 2015 European Synchrotron Radiation Facility, Grenoble, France
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
from ..id21_ffxas import process
from ...io.edf import saveedf
from ...common.listtools import move

from testfixtures import TempDirectory

import os

import numpy as np
import h5py

class test_alignSource(unittest.TestCase):

    def setUp(self):
        self.dir = TempDirectory()

    def tearDown(self):
        self.dir.cleanup()

    def getdata(self,path,radix):
        n1,n2 = 55,65

        xv = range(int(n2/3),int(2*n2/3))
        yv = range(int(2*n2/3),int(n2/3),-1)
        n = min(len(xv),3)
        energy = np.linspace(7,7.3,n) # keV
        intensity = np.linspace(200,100,n) # DU/sec
        transmission = np.linspace(0.9,0.3,n)
        nbdata = 8
        nbflat = 2
        tdata = 0.4 # sec
        tflat = 0.1 # sec
        
        darkoff = 50 # DU
        darkgain = 33 # DU/sec
        nbdark = 10
        
        darkdata = np.full((n1,n2),(darkoff+darkgain*tdata)*nbdark)
        hdark = {"energy":energy.max(),"exposure_time":tdata,"nb_frames":nbdark}
        filename = os.path.join(path,"{}_dark_{}_0000.edf".format(radix,tdata))
        saveedf(filename,darkdata,hdark)
        
        if tflat!=tdata:
            darkflat = np.full((n1,n2),(darkoff+darkgain*tflat)*nbdark)
            hdark = {"energy":energy.max(),"exposure_time":tflat,"nb_frames":nbdark}
            filename = os.path.join(path,"{}_dark_{}_0000.edf".format(radix,tflat))
            saveedf(filename,darkflat,hdark)
        
        for c,(x,y,e,i,t) in enumerate(zip(xv,yv,energy,intensity,transmission)):
            data = np.full((n1,n2),(darkoff + (i + darkgain)*tdata)*nbdata)
            data[y,x] = (darkoff + (t*i + darkgain)*tdata)*nbdata
            
            hdata = {"energy":e,"exposure_time":tdata,"nb_frames":nbdata}
            filename = os.path.join(path,"{}_data_0000_{:04d}_0000.edf".format(radix,c))
            saveedf(filename,data,hdata)
            
            flat1 = np.full((n1,n2),(darkoff + (i*1.05 + darkgain)*tflat)*nbflat)
            flat2 = np.full((n1,n2),(darkoff + (i*0.95 + darkgain)*tflat)*nbflat)
            hflat = {"energy":e,"exposure_time":tflat,"nb_frames":nbflat}
            filename = os.path.join(path,"{}_ref_0000_{:04d}_{:04d}.edf".format(radix,c,0))
            saveedf(filename,flat1,hflat)
            filename = os.path.join(path,"{}_ref_0000_{:04d}_{:04d}.edf".format(radix,c,1))
            saveedf(filename,flat2,hflat)
        
        
        return {"energy":energy,\
                "intensity":intensity,\
                "transmission":transmission,\
                "xv":xv,\
                "yv":yv,\
                "n1":n1,\
                "n2":n2,\
                "n":n}
    
    def align(self,sourcepath,radix,outname,ext,alignmethod,refimageindex,roiraw,roialign,roiresult,crop,stackdim):
        # Raw data
        rebin = (1,1)
        skippreprocessing = False

        # Result
        destpath = os.path.join(sourcepath,"results",outname)

        # Normalization
        skipnormalization = False

        # Alignment
        skipalign = False
        plot = False

        # Process
        process(sourcepath,destpath,radix,ext,rebin,alignmethod,\
            skippre=skippreprocessing,skipnormalization=skipnormalization,skipalign=skipalign,\
            roiraw=roiraw,roialign=roialign,roiresult=roiresult,\
            refimageindex=refimageindex,crop=crop,plot=plot,stackdim=stackdim)
    
    def checkresult(self,sourcepath,outname,params):
        destpath = os.path.join(sourcepath,"results",outname)

        # Check normalized results
        imgshape = (params["n1"],params["n2"])
        shape = tuple(move([params["n1"],params["n2"],params["n"]],2,params["stackdim"]))
        
        with h5py.File(os.path.join(destpath,"{}.h5".format(outname))) as f:
            np.testing.assert_allclose(params["energy"],f["detector0"]["sample"]["energy"])
            np.testing.assert_array_equal(range(params["n1"]),f["detector0"]["sample"]["row"])
            np.testing.assert_array_equal(range(params["n2"]),f["detector0"]["sample"]["col"])
            
            fdata = f["detector0"]["sample"]["data"]
            self.assertEqual(shape,fdata.shape)
            for i,(t,x,y) in enumerate(zip(params["transmission"],params["xv"],params["yv"])):
                index = tuple(move([slice(None),slice(None),i],2,params["stackdim"]))
                data = f["detector0"]["sample"]["data"][index]
                
                data2 = np.zeros(imgshape,dtype=np.float32)
                data2[y,x] = -np.log(t)
                np.testing.assert_allclose(data,data2,atol=1e-6)

        # Check aligned results
        off = [1,-1]
        totaloff = [off[0]*(params["n"]-1),off[1]*(params["n"]-1)]

        if params["crop"]:
            imgshape = (params["n1"]-abs(totaloff[0]),params["n2"]-abs(totaloff[1]))
        else:
            imgshape = (params["n1"]+abs(totaloff[0]),params["n2"]+abs(totaloff[1]))
        shape = [imgshape[0],imgshape[1],params["n"]]
        shape = tuple(move(shape,2,params["stackdim"]))
        
        with h5py.File(os.path.join(destpath,"{}.align.h5".format(outname))) as f:
            np.testing.assert_allclose(params["energy"],f["detector0"]["sample"]["energy"])
            
            if params["crop"]:
                row = np.arange(max(totaloff[0],0),params["n1"]+min(totaloff[0],0))
                col = np.arange(max(totaloff[1],0),params["n2"]+min(totaloff[1],0))
            else:
                row = np.arange(min(totaloff[0],0),params["n1"]+max(totaloff[0],0))
                col = np.arange(min(totaloff[1],0),params["n2"]+max(totaloff[1],0))
            np.testing.assert_array_equal(row,f["detector0"]["sample"]["row"])
            np.testing.assert_array_equal(col,f["detector0"]["sample"]["col"])
            
            fdata = f["detector0"]["sample"]["data"]
            self.assertEqual(shape,fdata.shape)
            
            y = np.nanargmin(np.abs(row-params["yv"][0]))
            x = np.nanargmin(np.abs(col-params["xv"][0]))
            data2 = np.zeros(imgshape,dtype=np.float32)
            
            for i,t in enumerate(params["transmission"]):
                index = tuple(move([slice(None),slice(None),i],2,params["stackdim"]))
                data = f["detector0"]["sample"]["data"][index]
                data[np.isnan(data)] = 0
                     
                data2[y,x] = -np.log(t)
                np.testing.assert_allclose(data,data2,atol=1e-6)
                
            trn = np.arange(params["n"])
            trn = np.stack([trn,-trn],axis=1)
            np.testing.assert_allclose(f["processing"]["2.align"]["changeofframes"][:],trn)
    
    def test_process(self):
        sourcepath = self.dir.path
        radix = "ff"
        params = self.getdata(sourcepath,radix)
        
        roiraw = None
        roialign = None
        roiresult = None
        ext = ""
        alignmethod = "max"
        refimageindex = None
        outname = radix
        roialign = ((3,-3),(4,-4))

        for crop in [True,False]:
            for roialign in [None,((3,-3),(4,-4))]:
                for stackdim in [0,1,2]:
                    self.align(sourcepath,radix,outname,ext,alignmethod,refimageindex,roiraw,roialign,roiresult,crop,stackdim)
                    params["stackdim"] = stackdim
                    params["crop"] = crop
                    self.checkresult(sourcepath,outname,params)

def test_suite_all():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_alignSource("test_process"))
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite_all()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
