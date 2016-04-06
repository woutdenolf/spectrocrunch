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
import sift_pyocl as sift
import os
import numpy as np

class alignSift(align):

    def __init__(self,*args,**kwargs):
        super(alignSift,self).__init__(*args,**kwargs)

        # pyopencl stuff
        device = sift.opencl.ocl.select_device(type="CPU", best=True)
        self.ctx = sift.opencl.pyopencl.Context(devices=[sift.opencl.pyopencl.get_platforms()[device[0]].get_devices()[device[1]]])
        self.queue = sift.opencl.pyopencl.CommandQueue(self.ctx)
        
        # Prepare alignment kernel
        self.max_workgroup_size = None
        sift.param.par.Scales = 8
        sift.param.par.PeakThresh = 0.
        self.inshape = self.source.imgsize
        self.outshape = self.source.imgsize

        self.newsiftplan()
        self.matchplan = sift.MatchPlan(context=self.ctx, max_workgroup_size=self.max_workgroup_size)
        self.kp1 = None
        self.kp2 = None
        
        # Prepare transformation kernel
        self.workgroupshape = (8, 4)
        kernel_file = os.path.join(os.path.dirname(os.path.abspath(sift.__file__)),"transform.cl")
        kernel_src = open(kernel_file).read()
        self.transformix = sift.opencl.pyopencl.Program(self.ctx, kernel_src).build()#('-D WORKGROUP_SIZE=%s' % self.max_workgroup_size)
        self.newtransformixshape()

        # Prepare transformation buffers
        self.buffers = {}
        self.newtransformationIObuffer()
        self.offset = self.idoffset.copy()
        self.linear = self.idlinear.copy()
        self.buffers["matrix"] = sift.opencl.pyopencl.array.empty(self.queue, shape=(2, 2), dtype=self.dtype)
        self.buffers["offset"] = sift.opencl.pyopencl.array.empty(self.queue, shape=(1, 2), dtype=self.dtype)
        cpy1 = sift.opencl.pyopencl.enqueue_copy(self.queue, self.buffers["matrix"].data, self.linear)
        cpy2 = sift.opencl.pyopencl.enqueue_copy(self.queue, self.buffers["offset"].data, self.offset)

    def newsiftplan(self):
        """New kernel for finding keypoints
        """
        self.siftplan = sift.SiftPlan(shape = self.inshape, dtype = self.dtype, context=self.ctx, max_workgroup_size=self.max_workgroup_size)

    def newtransformationIObuffer(self):
        """New IO buffers for the transformation kernel
        """
        self.buffers["input"] = sift.opencl.pyopencl.array.empty(self.queue, shape=self.inshape, dtype=self.dtype)
        self.buffers["output"] = sift.opencl.pyopencl.array.empty(self.queue, shape=self.outshape, dtype=self.dtype)

    def newtransformixshape(self):
        shape = self.inshape[::-1]
        self.transformixshape = tuple((int(i) + int(j) - 1) & ~(int(j) - 1) for i, j in zip(shape, self.workgroupshape))

    def changeshape(self,shape):
        """Adapt shape dependent buffers and kernels
        """
        if self.inshape == shape:
            return
        self.inshape = shape
        self.outshape = shape
        self.newtransformationIObuffer()
        self.newtransformixshape()

    def changerefshape(self,shape):
        """Adapt shape dependent buffers and kernels
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
        if self.isidentity():
            return img
        self.changeshape(img.shape)

        # Copy image to buffer
        data = np.ascontiguousarray(img, self.dtype)
        cpy = sift.opencl.pyopencl.enqueue_copy(self.queue, self.buffers["input"].data, data)
        cpy.wait()

        # Apply transformation
        self.execute_transformatrix()
        return self.buffers["output"].get()
        
    def execute_alignkernel(self,img):
        """Align image on reference
        """

        # Copy image to buffer
        data = np.ascontiguousarray(img, self.dtype)
        cpy = sift.opencl.pyopencl.enqueue_copy(self.queue, self.buffers["input"].data, data)
        cpy.wait()

        # Find keypoints of buffered image
        self.kp2 = self.siftplan.keypoints(self.buffers["input"])

        # Find matching reference keypoints
        self.offset[:] = 0
        if self.kp1.size != 0 and self.kp2.size != 0:
            raw_matching = self.matchplan.match(self.buffers["ref_kp_gpu"], self.kp2, raw_results=True)

            # Extract transformation from matching keypoints
            matching = np.recarray(shape=raw_matching.shape, dtype=self.matchplan.dtype_kp)
            len_match = raw_matching.shape[0]
            if len_match != 0:
                matching[:, 0] = self.kp1[raw_matching[:, 0]]
                matching[:, 1] = self.kp2[raw_matching[:, 1]]

                # Extract transformation matrix
                dx = matching[:, 1].x - matching[:, 0].x
                dy = matching[:, 1].y - matching[:, 0].y
                self.offset[:] = (np.median(dy),np.median(dx)) # y is the first dimension in python
                

        # Apply transformation
        if self.offset[0]==0 and self.offset[1]==0:
            return img

        #cpy1 = sift.opencl.pyopencl.enqueue_copy(self.queue, self.buffers["matrix"].data, self.linear)
        cpy2 = sift.opencl.pyopencl.enqueue_copy(self.queue, self.buffers["offset"].data, self.offset)
        self.execute_transformatrix()
        return self.buffers["output"].get()

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

    def isidentity(self):
        """Transformation is the identity?
        """
        return np.array_equal(self.offset,self.idoffset) and np.array_equal(self.linear,self.idlinear)

    def set_reference(self,img,previous=False):
        """Reference for alignment
        """
        self.changerefshape(img.shape)
        
        # Set reference keypoints
        if previous:
            self.kp1 = self.kp2
        else:
            self.kp1 = self.siftplan.keypoints(img)
        self.buffers["ref_kp_gpu"] = sift.opencl.pyopencl.array.to_device(self.matchplan.queue, self.kp1)

    def get_transformation(self):
        """Get transformation from alignment kernel
        """
        return self.offset

    def set_transformation(self,offset,bchanged):
        """Set the transformation kernel according to the alignment kernel and adapted transformation
        """
        if bchanged:
            self.offset[:] = offset
            #cpy1 = sift.opencl.pyopencl.enqueue_copy(self.queue, self.buffers["matrix"].data, self.linear)
            cpy2 = sift.opencl.pyopencl.enqueue_copy(self.queue, self.buffers["offset"].data, self.offset)
