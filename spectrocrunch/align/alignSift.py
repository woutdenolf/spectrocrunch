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


from .align import align

from silx.image import sift
from silx.opencl import ocl
from silx.opencl.utils import get_opencl_code
import pyopencl

import os
import numpy as np
from scipy import stats
from .types import transformationType
import logging

class alignSift(align):

    def __init__(self,*args,**kwargs):
        super(alignSift,self).__init__(*args,**kwargs)

        # No 1D
        if 1 in self.source.imgsize:
            raise ValueError("Sift can only be applied on images, not 1D vectors.")

        # pyopencl stuff
        device = ocl.select_device(type="GPU", best=True)
        if device is None:
            device = ocl.select_device(best=True)
        self.ctx = pyopencl.Context(devices=[pyopencl.get_platforms()[device[0]].get_devices()[device[1]]])
        self.queue = pyopencl.CommandQueue(self.ctx)
        
        # Prepare alignment kernel
        self.max_workgroup_size = None
        sift.param.par.Scales = 8
        sift.param.par.PeakThresh = 0.
        self.inshape = self.source.imgsize
        self.outshape = self.source.imgsize

        if min(self.inshape)<=5:
            sift.param.par["BorderDist"] = 0
        if 3*(2 * sift.param.par["BorderDist"] + 2) > min(self.inshape):
            sift.param.par["BorderDist"] = 1

        self.newsiftplan()
        self.matchplan = sift.MatchPlan(context=self.ctx, max_workgroup_size=self.max_workgroup_size)
        self.kp1 = None
        self.kp2 = None
        
        # Prepare transformation kernel
        self.workgroupshape = (8, 4)
        self.transformix = pyopencl.Program(self.ctx, get_opencl_code("transform.cl")).build()#('-D WORKGROUP_SIZE=%s' % self.max_workgroup_size)
        self.newtransformixshape()
        
        # Prepare transformation buffers
        self.buffers = {}
        self.newtransformationIObuffer()
        self._transform = self.defaulttransform()
        self.buffers["matrix"] = pyopencl.array.empty(self.queue, shape=(2, 2), dtype=self.dtype)
        self.buffers["offset"] = pyopencl.array.empty(self.queue, shape=(1, 2), dtype=self.dtype)
        self.updatecofbuffer()

    def updatecofbuffer(self):
        # (x,y) becomes (y,x)
        cpy1 = pyopencl.enqueue_copy(self.queue, self.buffers["matrix"].data, np.ascontiguousarray(self._transform.getlinear().T))
        cpy2 = pyopencl.enqueue_copy(self.queue, self.buffers["offset"].data, np.ascontiguousarray(self._transform.gettranslation()[::-1]))

    def newsiftplan(self):
        """New kernel for finding keypoints
        """
        self.siftplan = sift.SiftPlan(shape = self.inshape, dtype = self.dtype, context=self.ctx, max_workgroup_size=self.max_workgroup_size)

    def newtransformationIObuffer(self):
        """New IO buffers for the transformation kernel
        """
        self.buffers["input"] = pyopencl.array.empty(self.queue, shape=self.inshape, dtype=self.dtype)
        self.buffers["output"] = pyopencl.array.empty(self.queue, shape=self.outshape, dtype=self.dtype)

    def newtransformixshape(self):
        shape = self.inshape[::-1]
        self.transformixshape = tuple((int(i) + int(j) - 1) & ~(int(j) - 1) for i, j in zip(shape, self.workgroupshape))

    def changeshape(self,shape):
        """Adapt shape dependent buffers and kernels for transformation
        """
        if self.inshape == shape:
            return
        self.inshape = shape
        self.outshape = shape
        self.newtransformationIObuffer()
        self.newtransformixshape()

    def changerefshape(self,shape):
        """Adapt shape dependent buffers and kernels for alignment
        """
        if self.inshape == shape:
            return
        self.inshape = shape
        self.outshape = shape
        self.newsiftplan()
        self.newtransformationIObuffer()

    def execute_transformkernel(self,img):
        """Transform image according with the transformation kernel
        """
        if self._transform.isidentity():
            return img
        self.changeshape(img.shape)

        # Copy image to buffer
        data = np.ascontiguousarray(img, self.dtype)
        cpy = pyopencl.enqueue_copy(self.queue, self.buffers["input"].data, data)
        cpy.wait()

        # Apply transformation
        self.execute_transformatrix()
        return self.buffers["output"].get()
        
    def execute_alignkernel(self,img):
        """Align image on reference
        """

        # Copy image to buffer
        data = np.ascontiguousarray(img, self.dtype)
        cpy = pyopencl.enqueue_copy(self.queue, self.buffers["input"].data, data)
        cpy.wait()

        # Find keypoints of buffered image
        self.kp2 = self.siftplan.keypoints(self.buffers["input"])

        # Find matching reference keypoints
        self._transform.settranslation(np.zeros(2))
        if self.kp1.size != 0 and self.kp2.size != 0:
            raw_matching = self.matchplan.match(self.buffers["ref_kp_gpu"], self.kp2, raw_results=True)

            # Extract transformation from matching keypoints
            matching = np.recarray(shape=raw_matching.shape, dtype=self.matchplan.dtype_kp)
            len_match = raw_matching.shape[0]
            if len_match != 0:
                matching[:, 0] = self.kp1[raw_matching[:, 0]]
                matching[:, 1] = self.kp2[raw_matching[:, 1]]

                # Transformation from matching keypoints
                self.transformationFromKp(matching[:, 0].x,matching[:, 0].y,matching[:, 1].x,matching[:, 1].y)

                # Extract transformation matrix
                #dx = matching[:, 1].x - matching[:, 0].x
                #dy = matching[:, 1].y - matching[:, 0].y
                #self._transform.settranslation((np.median(dy),np.median(dx))) # y is the first dimension in python
                
        # Apply transformation
        if self._transform.isidentity():
            return img

        self.updatecofbuffer()
        self.execute_transformatrix()
        return self.buffers["output"].get()

    def centroid(self,x):
        #return np.mean(x)
        #return np.median(x)
        return stats.trim_mean(x,0.1) # trimmed mean (trim 10% at both sides)

    def solvelinearsystem(self,A,b):
        ##### Using pseudo-inverse #####
        # A-1* = (A^T.A)^(-1).A^T
        try:
            S = np.dot(A.T, A)
            sol = np.dot(np.linalg.inv(S), np.dot(A.T, b))
        except np.linalg.LinAlgError as err:
            logger = logging.getLogger(__name__)
            logger.error("Singular matrix in calculating a transformation from SIFT keypoints")
            sol = None
            
        #sol = np.dot(numpy.linalg.pinv(A),b) #slower?

        ##### Using SVD #####
        #sol = np.linalg.lstsq(A,b)[0] # computing the numpy solution

        ##### Using QR #####
        #Q,R = np.linalg.qr(A) # qr decomposition of A
        #Qb = np.dot(Q.T,b) # computing Q^T*b (project b onto the range of A)
        #sol = np.linalg.solve(R,Qb) # solving R*x = Q^T*b

        # result
        #MSE = np.linalg.norm(b - np.dot(A,sol))**2/N #Mean Squared Error
        return sol

    def transformationFromKp(self,xsrc,ysrc,xdest,ydest):
        """ Least-squares transformation parameters to map src to dest

            Remark: the rigid transformation is the most problematic (cfr. test_sift_mapping)
        """
        self._transform.setidentity()

        if self.transfotype==transformationType.translation:
            self._transform.settranslation([self.centroid(xdest-xsrc),self.centroid(ydest-ysrc)])

        elif self.transfotype==transformationType.rigid:
            # https://igl.ethz.ch/projects/ARAP/svd_rot.pdf
            #
            # R.(X-Xcen) = Y-Ycen
            # R.X + T = Y   and   T = Ycen - R.Xcen

            censrc = np.asarray([self.centroid(xsrc),self.centroid(ysrc)])
            cendest = np.asarray([self.centroid(xdest),self.centroid(ydest)])

            XT = np.column_stack((xsrc,ysrc)) 
            YT = np.column_stack((xdest,ydest))

            S = np.dot(np.dot(np.transpose(XT-censrc),np.identity(len(xsrc))),YT-cendest)
            U, s, V = np.linalg.svd(S, full_matrices=False)
            C = np.dot(V,np.transpose(U))
            self._transform.setlinear(C)

            YTnoshift = YT-np.dot(XT,C.T)
            #self._transform.settranslation(cendest - np.dot(C,censrc))
            # This seems to be more accurate
            self._transform.settranslation([self.centroid(YTnoshift[:,0]),self.centroid(YTnoshift[:,1])])
            
        elif self.transfotype==transformationType.similarity:
            # Similarity transformation:
            #    x' = a.x - b.y + t0
            #    y' = b.x + a.y + t1
            #    sol = [a,b,t0,t1]
            #
            # xsrc = [x1,x2,...]
            # ysrc = [y1,y2,...]
            # xdest = [x1',x2',...]
            # ydest = [y1',y2',...]
            #
            # X = x1 -y1  1  0
            #     y1  x1  0  1
            #     x2 -y2  1  0
            #     y2  x2  0  1
            #     ...
            #
            # Y = x1'
            #     y1'
            #     x2'
            #     y2'
            #     ...

            N = len(xsrc)
            
            X = np.zeros((2 * N, 4))
            X[::2, 0] = xsrc
            X[1::2, 0] = ysrc
            X[::2, 1] = -ysrc
            X[1::2, 1] = xsrc
            X[::2, 2] = 1
            X[1::2, 3] = 1

            Y = np.zeros((2 * N, 1))
            Y[::2, 0] = xdest
            Y[1::2, 0] = ydest

            sol = self.solvelinearsystem(X,Y)
            if sol is not None:
                self._transform.setlinear([[sol[0],-sol[1]],[sol[1],sol[0]]])
                self._transform.settranslation(sol[2:].flatten())

        elif self.transfotype==transformationType.affine:
            # Affine transformation:
            #    x' = a.x + b.y + t0
            #    y' = c.x + d.y + t1
            #    sol = [a,b,t0,c,d,t1]
            #
            # xsrc = [x1,x2,...]
            # ysrc = [y1,y2,...]
            # xdest = [x1',x2',...]
            # ydest = [y1',y2',...]
            #
            # X = x1 y1  1  0  0  0
            #      0  0  0 x1 y1  1
            #     x2 y2  1  0  0  0
            #      0  0  0 x2 y2  1
            #     ...
            #
            # Y = x1'
            #     y1'
            #     x2'
            #     y2'
            #     ...

            raise NotImplementedError("Sift doesn't support affine transformations (keypoints not invariant).")

            N = len(xsrc)

            X = np.zeros((2 * N, 6))
            X[::2, 0] = xsrc
            X[::2, 1] = ysrc
            X[::2, 2] = 1
            X[1::2, 3] = xsrc
            X[1::2, 4] = ysrc
            X[1::2, 5] = 1
            
            Y = np.zeros((2 * N, 1))
            Y[::2, 0] = xdest
            Y[1::2, 0] = ydest

            sol = self.solvelinearsystem(X,Y)
            if sol is not None:
                self._transform.setaffine(sol.reshape((2,3)))

        elif self.transfotype==transformationType.homography:
            # Projective transformation:
            #    x' = (a.x + b.y + t0)/(px.x+py.y+1)
            #    y' = (c.x + d.y + t1)/(px.x+py.y+1)
            #     x' = a.x + b.y + t0 - px.x.x' - py.y.x'
            #     y' = c.x + d.y + t1 - px.x.y' - py.y.y'
            #    sol = [a,b,t0,c,d,t1,px,py]
            #
            # xsrc = [x1,x2,...]
            # ysrc = [y1,y2,...]
            # xdest = [x1',x2',...]
            # ydest = [y1',y2',...]
            #
            # X = x1 y1  1  0  0  0 -x1.x1' -y1.x1' -x1'
            #      0  0  0 x1 y1  1 -x1.y1' -y1.y1' -y1'
            #     x2 y2  1  0  0  0 -x2.x2' -y2.x2' -x1'
            #      0  0  0 x2 y2  1 -x2.y2' -y2.y2' -y1'
            #     ...
            #
            # Y = 0
            #     0
            #     0
            #     0
            #     ...
            
            raise NotImplementedError("Sift doesn't support homographies (keypoints not invariant).")
            #if self.usekernel:
            #    raise NotImplementedError("Sift doesn't support this type of transformation.")

            N = len(xsrc)

            X = np.zeros((2 * N, 8))
            X[::2, 0] = xsrc
            X[::2, 1] = ysrc
            X[::2, 2] = 1
            X[1::2, 3] = xsrc
            X[1::2, 4] = ysrc
            X[1::2, 5] = 1
            X[::2, 6] = -xsrc*xdest
            X[1::2, 6] = -xsrc*ydest
            X[::2, 7] = -ysrc*xdest
            X[1::2, 7] = -ysrc*ydest

            Y = np.zeros((2 * N, 1))
            Y[::2, 0] = xdest
            Y[1::2, 0] = ydest

            sol = self.solvelinearsystem(X,Y)
            if sol is not None:
                self._transform.sethomography(np.append(sol,1).reshape((3,3))) 
        else:
            raise NotImplementedError("Sift doesn't support this type of transformation.")

    def execute_transformatrix(self):
        """Execute transformation kernel
        """
        ev = self.transformix.transform(self.queue, self.transformixshape, self.workgroupshape,
                                   self.buffers["input"].data,
                                   self.buffers["output"].data,
                                   self.buffers["matrix"].data,
                                   self.buffers["offset"].data,
                                   np.int32(self.inshape[1]), # input width
                                   np.int32(self.inshape[0]), # input height
                                   np.int32(self.outshape[1]), # output width
                                   np.int32(self.outshape[0]), # output height
                                   np.float32(self.cval), # fill
                                   #self.siftplan.buffers["min"].get()[0],
                                   np.int32(1)) # bilinear interpolation

    def set_reference(self,img,previous=False):
        """Reference for alignment
        """
        self.changerefshape(img.shape)
        
        # Set reference keypoints
        if previous:
            self.kp1 = self.kp2
        else:
            self.kp1 = self.siftplan.keypoints(np.ascontiguousarray(img, self.dtype))

        self.buffers["ref_kp_gpu"] = pyopencl.array.to_device(self.matchplan.queue, self.kp1)

    def get_transformation(self):
        """Get transformation from alignment kernel
        """
        return self._transform

    def set_transformation(self,transform,bchanged):
        """Set the transformation kernel according to the alignment kernel and adapted transformation
        """
        if bchanged:
            self._transform.set(transform)
            self.updatecofbuffer()

