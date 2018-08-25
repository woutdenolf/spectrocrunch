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

from ..transform import transform
from ..types import transformationType
from .import helper_teststack

import numpy as np
from skimage import data
import scipy.ndimage.interpolation as ndtransform
import skimage.transform as sktransform

class test_transform(unittest.TestCase):

    def genkeypoints(self,N=10):
        N = 10
        xsrc = np.random.random(N)*100
        ysrc = np.random.random(N)*100
        return np.vstack((xsrc,ysrc,np.ones(N)))
    
    def getxy(self,Y):
        return Y[0,:]/Y[2,:],Y[1,:]/Y[2,:]
    
    def test_linearmapping(self):

        # Generate points
        X = self.genkeypoints(N=10)
        xsrc,ysrc = self.getxy(X)
        
        # TODO: improve on rigid transformation
        types = [transformationType.translation, transformationType.rigid, transformationType.similarity, transformationType.affine]
        rtols = [1e-7, 0.5, 1e-5, 1e-5]
        atols = [1e-6, 0.5, 1e-5, 1e-5]

        for i in range(10):
            for t,rtol,atol in zip(types,rtols,atols):
                Mcoord,Mcof = helper_teststack.transformation(t,2)

                # Map src to dest
                Y = Mcoord.dot(X)
                xdest,ydest = self.getxy(Y)

                # Add noise
                #xdest += np.random.random(N)-0.5
                #ydest += np.random.random(N)-0.5

                # Get transformation
                o = transform(t,dtype=X.dtype,cval=np.nan)
                o.fromkeypoints(xsrc,ysrc,xdest,ydest)
                Mcof2 = o.getnumpyhomography()
                
                np.testing.assert_allclose(Mcof,Mcof2,rtol=0,atol=atol)
                np.testing.assert_allclose(X,Mcof2.dot(Y),rtol=rtol,atol=0)
                  
    def test_compare(self):
        
        # REMARK: sktransform.warp needs change-of-frame (i.e. map dest keypoints -> src keypoints)
        
        o1 = sktransform.EuclideanTransform(translation=[10,20],rotation=np.radians(0.1))
        
        o2 = sktransform.EuclideanTransform()
        src = np.vstack([[0,0.],[5,3.]])
        dest = src + np.array([[10,20]])
        o2.estimate(dest,src)
        
        kwargs = {"cval":0,"order":1,"mode":"constant"}
        img = data.camera()[::-1,:]
        
        for o in [o1]:
            for i in range(100):
                cof = o.params
                # shift: takes transformation vector for coordinates (y,x)
                img1 = ndtransform.shift(img,-cof[1::-1,2],**kwargs)
                # affine_transform: takes change-of-frame matrix for coordinates (y,x)
                img2 = ndtransform.affine_transform(img,cof[0:2,0:2].T,offset=cof[1::-1,2],**kwargs)
                # warp: takes change-of-frame matrix
                img3 = sktransform.warp(img,o,**kwargs)

            if np.isclose(o.rotation,0):
                np.allclose(img1,img2)
            np.allclose(img2,img3)

            if False:
                import matplotlib.pyplot as plt
                fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, figsize=(8, 4),
                                               sharex=True, sharey=True)
                ax1.imshow(img1,origin="lower")
                ax2.imshow(img2,origin="lower")
                ax3.imshow(img3,origin="lower")
                plt.show()
        
def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_transform("test_linearmapping")) 
    testSuite.addTest(test_transform("test_compare")) 
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)