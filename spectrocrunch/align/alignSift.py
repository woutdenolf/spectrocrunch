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
import warnings
import os
import numpy as np
import logging
try:
    import pyopencl
except ImportError:
    pyopencl = None
    warnings.warn("pyopencl is not installed", ImportWarning)

from .align import align
from .types import transformationType
from ..math.fit1d import lstsq

from silx.image import sift
from silx.opencl import ocl
from silx.opencl.utils import get_opencl_code


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
        self.siftdtype = np.float32
        sift.param.par.Scales = 8
        sift.param.par.PeakThresh = 0.
        self.inshape = self.source.imgsize
        self.outshape = self.source.imgsize

        if min(self.inshape)<=5:
            sift.param.par["BorderDist"] = 0
        if 3*(2 * sift.param.par["BorderDist"] + 2) > min(self.inshape):
            sift.param.par["BorderDist"] = 1

        self.newsiftplan()
        self.matchplan = sift.MatchPlan(ctx=self.ctx)
        self.kpref = None
        self.kpmoving = None
        
        # Prepare transformation kernel
        self.workgroupshape = (8, 4)
        try:
            self.transformix = pyopencl.Program(self.ctx, get_opencl_code("transform.cl")).build()
        except IOError:
            self.transformix = pyopencl.Program(self.ctx, get_opencl_code("sift/transform.cl")).build()
        self.newtransformixshape()
        
        # Prepare transformation buffers
        self.buffers = {}
        self.newtransformationIObuffer()
        self._transform = self.defaulttransform(dtype=self.siftdtype)
        self.buffers["matrix"] = pyopencl.array.empty(self.queue, shape=(2, 2), dtype=self.siftdtype)
        self.buffers["offset"] = pyopencl.array.empty(self.queue, shape=(1, 2), dtype=self.siftdtype)
        self.updatecofbuffer()

    def updatecofbuffer(self):
        # (x,y) becomes (y,x)
        cpy1 = pyopencl.enqueue_copy(self.queue, self.buffers["matrix"].data, np.ascontiguousarray(self._transform.getlinear().T))
        cpy2 = pyopencl.enqueue_copy(self.queue, self.buffers["offset"].data, np.ascontiguousarray(self._transform.gettranslation()[::-1]))

    def newsiftplan(self):
        """New kernel for finding keypoints
        """
        self.siftplan = sift.SiftPlan(shape = self.inshape, dtype = self.siftdtype, ctx=self.ctx)

    def newtransformationIObuffer(self):
        """New IO buffers for the transformation kernel
        """
        self.buffers["input"] = pyopencl.array.empty(self.queue, shape=self.inshape, dtype=self.siftdtype)
        self.buffers["output"] = pyopencl.array.empty(self.queue, shape=self.outshape, dtype=self.siftdtype)

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
        data = np.ascontiguousarray(img, self.siftdtype)
        cpy = pyopencl.enqueue_copy(self.queue, self.buffers["input"].data, data)
        cpy.wait()

        # Apply transformation
        self.execute_transformatrix()
        return self.buffers["output"].get()
        
    def execute_alignkernel(self,img):
        """Align image on reference
        """

        # Copy image to buffer
        data = np.ascontiguousarray(img, self.siftdtype)
        cpy = pyopencl.enqueue_copy(self.queue, self.buffers["input"].data, data)
        cpy.wait()

        # Find keypoints of buffered image
        self.kpmoving = self.siftplan.keypoints(self.buffers["input"])

        # Find matching reference keypoints
        self._transform.settranslation(np.zeros(2))
        if self.kpref.size != 0 and self.kpmoving.size != 0:
            raw_matching = self.matchplan.match(self.buffers["ref_kp_gpu"], self.kpmoving, raw_results=True)

            # Extract transformation from matching keypoints
            matching = np.recarray(shape=raw_matching.shape, dtype=self.matchplan.dtype_kp)
            len_match = raw_matching.shape[0]
            if len_match != 0:
                matching[:, 0] = self.kpmoving[raw_matching[:, 1]]
                matching[:, 1] = self.kpref[raw_matching[:, 0]]
 
                # Map kpmoving to kpref
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

    def transformationFromKp(self,xsrc,ysrc,xdest,ydest):
        """ Least-squares transformation parameters to map src to dest

            Remark: the rigid transformation is the most problematic (cfr. test_sift_mapping)
        """
        self._transform.fromkeypoints(xsrc,ysrc,xdest,ydest)

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
                                   self.siftdtype(self.cval), # fill
                                   #self.siftplan.buffers["min"].get()[0],
                                   np.int32(1)) # bilinear interpolation

    def set_reference(self,img,previous=False):
        """Reference for alignment
        """
        self.changerefshape(img.shape)
        
        # Set reference keypoints
        if previous:
            self.kpref = self.kpmoving
        else:
            self.kpref = self.siftplan.keypoints(np.ascontiguousarray(img, self.siftdtype))

        self.buffers["ref_kp_gpu"] = pyopencl.array.to_device(self.matchplan.queue, self.kpref)

    def get_transformation(self):
        """Get transformation from alignment kernel
        """
        return self._transform

    def set_transformation(self,transform,bchanged):
        """Set the transformation kernel according to the alignment kernel and adapted transformation
        """
        if bchanged:
            self._transform.fromtransform(transform)
            self.updatecofbuffer()

