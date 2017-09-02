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

from testfixtures import TempDirectory

import os

import numpy as np

class test_alignSource(unittest.TestCase):

    def setUp(self):
        self.dir = TempDirectory()

    def tearDown(self):
        self.dir.cleanup()
        
    def getdata(self,path,radix):
        n1,n2 = 55,65

        xv = range(int(n2/3),int(2*n2/3))
        yv = range(int(n2/3),int(2*n2/3))
        energy = np.linspace(7,7.3,len(xv))
        intensity = np.linspace(200,100,len(xv)) # cts/sec
        att = np.linspace(1,0.3,len(xv))
        nbdata = 8
        nbflat = 2
        tdata = 0.4 # sec
        tflat = 0.1
        
        darkoff = 50
        darkgain = 33
        ndark = 10
        
        darkdata = np.full((n1,n2),(darkoff+darkgain*tdata)*ndark)
        hdark = {"energy":energy.max(),"exposure_time":tdata,"nb_frames":ndark}
        filename = os.path.join(path,"{}_dark_{}_0000.edf".format(radix,tdata))
        saveedf(filename,darkdata,hdark)
        
        flatm = 10/11.
        
        if tflat!=tdata:
            darkflat = np.full((n1,n2),(darkoff+darkgain*tflat)*ndark)
            hdark = {"energy":energy.max(),"exposure_time":tflat,"nb_frames":ndark}
            filename = os.path.join(path,"{}_dark_{}_0000.edf".format(radix,tdata))
            saveedf(filename,darkdata,hdark)
        
        for c,(x,y,e,i,a) in enumerate(zip(xv,yv,energy,intensity,att)):
            data = np.zeros((n1,n2))
            data[y,x] = (darkoff + (a*i + darkgain)*tdata)*nbdata
            hdata = {"energy":e,"exposure_time":tdata,"nb_frames":nbdata}
            filename = os.path.join(path,"{}_data_0000_{:04d}_0000.edf".format(radix,c))
            saveedf(filename,data,hdata)
            
            flat1 = np.full((n1,n2),(darkoff + (i*flatm + darkgain)*tflat)*nbflat)
            flat2 = np.full((n1,n2),(darkoff + (i/flatm + darkgain)*tflat)*nbflat)
            hflat = {"energy":e,"exposure_time":tflat,"nb_frames":nbflat}
            filename = os.path.join(path,"{}_ref_0000_{:04d}_{:04d}.edf".format(radix,c,0))
            saveedf(filename,flat1,hflat)
            filename = os.path.join(path,"{}_ref_0000_{:04d}_{:04d}.edf".format(radix,c,1))
            saveedf(filename,flat2,hflat)
    
    def align(self,sourcepath,radix,outname,ext,alignmethod,refimageindex,roiraw,roialign,roiresult):
        path = sourcepath

        # Raw data
        rebin = (1,1)
        skippreprocessing = False

        # Result
        destpath = os.path.join(path,"results",outname)

        # Normalization
        skipnormalization = False

        # Alignment
        skipalign = False
        crop = True
        plot = True

        # Process
        process(sourcepath,destpath,radix,ext,rebin,alignmethod,\
            skippre=skippreprocessing,skipnormalization=skipnormalization,skipalign=skipalign,\
            roiraw=roiraw,roialign=roialign,roiresult=roiresult,\
            refimageindex=refimageindex,crop=crop,plot=plot)
    
    def test_process(self):
        sourcepath = self.dir.path
        radix = "ff"
        self.getdata(sourcepath,radix)
        roiraw = None
        roialign = None
        roiresult = None
        ext = ""
        alignmethod = "sift"
        refimageindex = None
        outname = radix
        roialign = ((1,-1),(2,-2))
        self.align(sourcepath,radix,outname,ext,alignmethod,refimageindex,roiraw,roialign,roiresult)

def test_suite_all():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    #testSuite.addTest(test_alignSource("test_process"))
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite_all()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
