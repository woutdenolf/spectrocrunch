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

# OpenCL Execution Model:
#  work item: basic unit of work
#  kernel: the code that runs on a work item (think C function)
#  program: collection of kernels (think dynamic shared library)
#  command queue: execution queue of kernels and memory operations
#
# Problem space in data parallelism:
#  global shape: nD-grid of work items
#  local shape: subspace of work items that shares local memory
#               and can communicate
#
# Enqueue kernel:
#  local shape: limited by device (due to hardware units)
#               limited by kernel (due to memory usage)
#  global shape: must contain problem space
#                must be a multiple of local shape


def local_shape(kernel, device, kernellimit=None, msplit=64):
    # Device limit: for example
    #  WORK_GROUP_SIZE = 1024
    #  max_work_item_sizes = (64,64,64)
    #  -> (64,8,2)
    query = pyopencl.kernel_work_group_info.WORK_GROUP_SIZE
    wgsize = kernel.get_work_group_info(query, device)
    wgshapemax = device.max_work_item_sizes
    ndmax = len(wgshapemax)
    # Kernel limit (use LOCAL_MEM_SIZE?):
    if kernellimit:
        wgsize = min(wgsize, kernellimit)
    # Size to shape
    if ndmax == 1:
        msplit = int(min(msplit, wgshapemax[0]))
    else:
        msplit = int(min(msplit, *wgshapemax[:2]))
    if wgsize > msplit:
        wgshape = (msplit, wgsize // msplit)
    elif wgsize > 32:
        wgshape = (32, wgsize // 32)
    elif wgsize > 4:
        wgshape = (wgsize // 4, 4)
    else:
        wgshape = (wgsize, 1)
    # Device limit:
    return tuple(min(n, m) for n, m in zip(wgshape, wgshapemax))


def global_shape(problemshape, localshape):
    # Smallest multiple of local shape that contains problem space:
    return tuple((int(nmin) + int(nsub) - 1) & ~(int(nsub) - 1)
                 for nmin, nsub in zip(problemshape, localshape))


def create_context():
    ctx = None
    if ocl is None:
        return ctx
    platdev = ocl.select_device(type="GPU", best=True)
    if platdev is None:
        platdev = ocl.select_device(best=True)
    if platdev is not None:
        try:
            device = pyopencl.get_platforms()[platdev[0]].get_devices()[platdev[1]]
            ctx = pyopencl.Context(devices=[device])
        except pyopencl.RuntimeError:
            pass
    if ctx is None:
        for platform in pyopencl.get_platforms():
            for device in platform.get_devices():
                try:
                    ctx = pyopencl.Context(devices=[device])
                except pyopencl.RuntimeError:
                    pass
    return ctx


class alignSift(align):

    def __init__(self,*args,**kwargs):
        super(alignSift,self).__init__(*args,**kwargs)

        # No 1D
        if 1 in self.source.imgsize:
            raise ValueError("Sift can only be applied on images, not 1D vectors.")

        # pyopencl stuff
        self.ctx = create_context()
        self.queue = pyopencl.CommandQueue(self.ctx)
        
        # Prepare alignment kernel
        self.siftdtype = np.float32
        sift.param.par.Scales = 8
        sift.param.par.PeakThresh = 0.
        self.inshape = None
        self.outshape = None
        self.buffers = {}
        self.matchplan = sift.MatchPlan(ctx=self.ctx)
        self.kpref = None
        self.kpmoving = None
        
        # Prepare transformation kernel
        try:
            self.transformix = pyopencl.Program(self.ctx,
                                get_opencl_code("transform.cl")).build()
        except IOError:
            self.transformix = pyopencl.Program(self.ctx,
                                get_opencl_code("sift/transform.cl")).build()

        # Prepare transformation buffers
        self._transform = self.defaulttransform(dtype=self.siftdtype)
        self.buffers["matrix"] = pyopencl.array.empty(self.queue,
                                    shape=(2, 2), dtype=self.siftdtype)
        self.buffers["offset"] = pyopencl.array.empty(self.queue,
                                    shape=(1, 2), dtype=self.siftdtype)
        self.updatecofbuffer()

        # Problem space
        self.problemshape = self.source.imgsize

    @property
    def device(self):
        return self.ctx.devices[0]

    def updatecofbuffer(self):
        # (x,y) becomes (y,x)
        cpy1 = pyopencl.enqueue_copy(self.queue, self.buffers["matrix"].data, np.ascontiguousarray(self._transform.getlinear().T))
        cpy2 = pyopencl.enqueue_copy(self.queue, self.buffers["offset"].data, np.ascontiguousarray(self._transform.gettranslation()[::-1]))

    @property
    def problemshape(self):
        return self.inshape
    
    @problemshape.setter
    def problemshape(self, value):
        if value == self.inshape:
            return
        self.inshape = value
        self.outshape = value

        self.siftplan = sift.SiftPlan(shape=self.inshape,
                        dtype=self.siftdtype, ctx=self.ctx)
        if min(self.inshape)<=5:
            sift.param.par["BorderDist"] = 0
        if 3*(2 * sift.param.par["BorderDist"] + 2) > min(self.inshape):
            sift.param.par["BorderDist"] = 1

        self.buffers["input"] = pyopencl.array.empty(self.queue,
                        shape=self.inshape, dtype=self.siftdtype)
        self.buffers["output"] = pyopencl.array.empty(self.queue,
                        shape=self.outshape, dtype=self.siftdtype)

        localshape = local_shape(self.transformix.transform, self.device, kernellimit=128)
        self.transformix_localshape = localshape
        self.transformix_globalshape = global_shape(self.inshape, localshape)

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
            raw_matching = self.matchplan.match(self.buffers["ref_kp_gpu"],
                            self.kpmoving, raw_results=True)

            # Extract transformation from matching keypoints
            matching = np.recarray(shape=raw_matching.shape,
                        dtype=self.matchplan.dtype_kp)
            len_match = raw_matching.shape[0]
            if len_match != 0:
                matching[:, 0] = self.kpmoving[raw_matching[:, 1]]
                matching[:, 1] = self.kpref[raw_matching[:, 0]]
 
                # Map kpmoving to kpref
                self.transformationFromKp(matching[:, 0].x, matching[:, 0].y,
                                          matching[:, 1].x, matching[:, 1].y)

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
        ev = self.transformix.transform(self.queue,
                                   self.transformix_globalshape,
                                   self.transformix_localshape,
                                   self.buffers["input"].data,
                                   self.buffers["output"].data,
                                   self.buffers["matrix"].data,
                                   self.buffers["offset"].data,
                                   np.int32(self.inshape[1]),  # input width
                                   np.int32(self.inshape[0]),  # input height
                                   np.int32(self.outshape[1]),  # output width
                                   np.int32(self.outshape[0]),  # output height
                                   self.siftdtype(self.cval),  # fill
                                   np.int32(1))  # bilinear interpolation

    def set_reference(self, img, previous=False):
        """Reference for alignment
        """
        self.problemshape = img.shape
        
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

