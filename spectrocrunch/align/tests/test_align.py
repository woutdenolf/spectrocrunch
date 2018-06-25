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

from ..alignElastix import alignElastix
from ..alignSift import alignSift
from ..alignFFT import alignFFT
from ..alignSimple import alignMin
from ..alignSimple import alignMax
from ..alignSimple import alignCentroid
from ..alignSimple import alignGaussMax
from ..types import transformationType

from .teststack import teststack
from .teststack import gettransformedimage
from .teststack import transformation as gentransform

import numpy as np

class test_align(unittest.TestCase):
    def testrelativecof(self,cofs,cofrel,msg=None):
        for i in range(1,cofs.shape[0]):
            cofrelcalc = np.dot(np.linalg.inv(cofs[i-1,...]),cofs[i,...])
            np.testing.assert_almost_equal(cofrelcalc,cofrel,decimal=2,err_msg=msg)

    def test_align(self,alignclass,transfotype,realistic=False,subpixel=True):

        if transfotype==transformationType.translation and\
            alignclass!=alignSift and alignclass!=alignElastix:
            lst = [True,False]
        else:
            lst = [False]
    
        for vector in lst:
            for transposed in lst:
                if transposed and not vector:
                    continue

                # Prepare dataIO
                inputstack,cofrel,stackdim = teststack(transfotype,vector=vector,transposed=transposed,realistic = realistic,subpixel=subpixel)
                outputstack = [np.zeros(1,dtype=np.float32)]*len(inputstack)

                # References
                refdatasetindex = 0
                refimageindex = 0#len(inputstack)//2
                    
                # Prepare alignment
                o = alignclass(inputstack,None,outputstack,None,None,stackdim=stackdim,overwrite=True,plot=False,transfotype=transfotype)

                # Check alignment
                if vector:
                    if transposed:
                        roi = ((1,-3),(0,1))
                    else:
                        roi = ((0,1),(1,-3))
                else:
                    roi = ((1,-3),(1,-3))


                for i in range(4):
                    pad = (i & 1)==1
                    crop = (i & 2)==2

                    # Fixed reference
                    msg = "Alignment: Pad = {}, Crop = {}, 1D = {}, transposed = {}, type = {}".format(pad,crop,vector,transposed,"fixed")
                    o.align(refdatasetindex,refimageindex=refimageindex,pad = pad,crop = crop,roi=roi)
                    self.testrelativecof(o.absolute_cofs(homography=True),cofrel,msg=msg)
 
                    # Pairwise: align on raw
                    msg = "Alignment: Pad = {}, Crop = {}, 1D = {}, transposed = {}, type = {}".format(pad,crop,vector,transposed,"pairwise/raw")
                    o.align(refdatasetindex,onraw = True,pad = pad,crop = crop,roi=roi)
                    self.testrelativecof(o.absolute_cofs(homography=True),cofrel,msg=msg)

                    # Pairwise: align on aligned
                    msg = "Alignment: Pad = {}, Crop = {}, 1D = {}, transposed = {}, type = {}".format(pad,crop,vector,transposed,"pairwise")
                    o.align(refdatasetindex,onraw = False,pad = pad,crop = crop,roi=roi)
                    self.testrelativecof(o.absolute_cofs(homography=True),cofrel,msg=msg)

    def test_sift_mapping(self):
        # Initialize alignSift (not important)
        inputstack = [np.zeros((2,2,2),dtype=np.float32)]*5
        outputstack = [np.zeros(1,dtype=np.float32)]*5
        o = alignSift(inputstack,None,outputstack,None,None,stackdim=2,overwrite=True)

        # Generate points
        N = 20
        xsrc = np.random.random(N)*100
        ysrc = np.random.random(N)*100
        XT = np.column_stack((xsrc,ysrc,np.ones(N)))

        types = [transformationType.translation, transformationType.rigid, transformationType.similarity]
        for t in types:
            M,_ = gentransform(t,2)

            # Transform points
            YT = np.dot(XT,M.transpose())
            xdest = YT[:,0]/YT[:,2]
            ydest = YT[:,1]/YT[:,2]

            # Add noise
            #xdest += np.random.random(N)-0.5
            #ydest += np.random.random(N)-0.5

            # Get transformation
            o.transfotype = t
            o.transformationFromKp(xsrc,ysrc,xdest,ydest)
            if t==transformationType.rigid:
                np.testing.assert_almost_equal(M,o._transform.getnumpyhomography(),decimal=1)
            else:
                np.testing.assert_allclose(M,o._transform.getnumpyhomography())

    def test_fft_internals(self):
         # Initialize alignFFT (not important)
        inputstack = [np.zeros((2,2,2),dtype=np.float32)]*5
        outputstack = [np.zeros(1,dtype=np.float32)]*5
        o = alignFFT(inputstack,None,outputstack,None,None,stackdim=2,overwrite=True,transfotype=transformationType.similarity)

        # Test Fourier related things
        img = np.abs(np.fft.fft2(np.arange(7*8).reshape((7,8))))
        img2 = np.fft.ifftshift(np.fft.fftshift(img))
        np.testing.assert_allclose(img,img2)

        # Test log-polar
        mx = 10
        my = 10
        nx = 501
        ny = 401
        dx = 2*mx/(nx-1.)
        dy = 2*my/(ny-1.)
        xv, yv = np.meshgrid(np.linspace(-mx,mx,nx),np.linspace(-my,my,ny))

        sx = 2.
        sy = 1.
        angle = 0.
        data = [(0.,0.,sx,sy,angle,1000.)]
        fixed = gettransformedimage(xv,yv,data,angle=True).reshape(xv.shape)

        angle = -2*np.pi/180
        scale = 1.3#0.9
        sx *= scale
        sy *= scale
        data = [(0.,0.,sx,sy,angle,1000.)]
        moving = gettransformedimage(xv,yv,data,angle=True).reshape(xv.shape)

        a = scale*np.cos(angle)
        b = scale*np.sin(angle)
        R = np.array([[a,-b,0],[b,a,0],[0,0,1]])

        T = np.array([[1,0,mx/dx],[0,1,my/dy],[0,0,1]])
        Tinv = np.array([[1,0,-mx/dx],[0,1,-my/dy],[0,0,1]])

        M = np.dot(T,np.dot(R,Tinv))

        o.set_reference(fixed)
        aligned = o.execute_alignkernel(moving)
        
        

        #TODO: doesn't pass!
        np.testing.assert_almost_equal(M,o._transform.getnumpyhomography(),decimal=1)

    def test_elastix(self):
        types = [transformationType.translation]
        for t in types:
            self.test_align(alignElastix,t)

    def test_sift(self):
        types = [transformationType.translation, transformationType.rigid, transformationType.similarity]
        types = [transformationType.translation]
        for t in types:
            self.test_align(alignSift,t)

    def test_fft(self):
        types = [transformationType.translation]
        for t in types:
            self.test_align(alignFFT,t)

    def test_min(self):
        self.test_align(alignMin,transformationType.translation,realistic=True,subpixel=False)

    def test_max(self):
        self.test_align(alignMax,transformationType.translation,subpixel=False)

    def test_centroid(self):
        self.test_align(alignCentroid,transformationType.translation,subpixel=False)

    def test_gaussmax(self):
        self.test_align(alignGaussMax,transformationType.translation,subpixel=False)

def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    ##testSuite.addTest(test_align("test_sift_mapping")) # not working for now
    ##testSuite.addTest(test_align("test_fft_internals")) # not working for now
    testSuite.addTest(test_align("test_min"))
    testSuite.addTest(test_align("test_max"))
    testSuite.addTest(test_align("test_centroid"))
    testSuite.addTest(test_align("test_gaussmax"))
    testSuite.addTest(test_align("test_fft"))
    testSuite.addTest(test_align("test_sift"))
    testSuite.addTest(test_align("test_elastix"))
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
