# -*- coding: utf-8 -*-

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
from silx.image import sift
from silx.opencl.sift.alignment import calc_size


class alignSift(align):
    def __init__(self, *args, **kwargs):
        super(alignSift, self).__init__(*args, **kwargs)

        # No 1D
        if 1 in self.source.imgsize:
            raise ValueError("Sift can only be applied on images, not 1D vectors.")

        # Prepare alignment kernel
        self.siftdtype = np.float32  # hardcoded in silx
        # sift.param.par.Scales = 8
        # sift.param.par.PeakThresh = 0.0
        if self.transfotype == transformationType.translation:
            self._shift_only = True
        elif self.transfotype == transformationType.affine:
            self._shift_only = False
        else:
            raise NotImplementedError(
                "Elastix doesn't support this type transformation."
            )

        # change of reference frame
        self._transform = self.defaulttransform(dtype=self.siftdtype)

    @property
    def _gpu_queue(self):
        return self.helper.queue

    @property
    def _gpu_buffers(self):
        return self.helper.cl_mem

    def _transform_to_buffer(self):
        """Copy transformation to GPU buffer"""
        # Transformation has (x,y) coordinates
        # The GPU buffer expects (y,x) coordinates
        pyopencl.enqueue_copy(
            self._gpu_queue,
            self._gpu_buffers["matrix"].data,
            np.ascontiguousarray(self._transform.getlinear().T),
        )
        pyopencl.enqueue_copy(
            self._gpu_queue,
            self._gpu_buffers["offset"].data,
            np.ascontiguousarray(self._transform.gettranslation()[::-1]),
        )

    def _buffer_to_transform(self):
        """Copy GPU buffer to transformation"""
        T = self._gpu_buffers["offset"].get()[0, ::-1]
        if self._shift_only:
            self._transform.settranslation(*T)
        else:
            M = np.identity(4)
            M[:2, 2] = T
            M[:2, :2] = self._gpu_buffers["matrix"].get().T
            self._transform.setaffine(M)

    def execute_transformkernel(self, img):
        """Transform image according with the transformation kernel
        """
        if self._transform.isidentity():
            return img

        # Copy image to GPU buffer
        data = np.ascontiguousarray(img, self.siftdtype)
        cpy = pyopencl.enqueue_copy(
            self._gpu_queue, self._gpu_buffers["input"].data, data
        )
        cpy.wait()

        # Apply transformation
        self.execute_transformix()
        return self._gpu_buffers["output"].get()

    def execute_alignkernel(self, img):
        """Align image on reference
        """
        result = self.helper.align(img, shift_only=self._shift_only)
        self._buffer_to_transform()
        return result

    def execute_transformix(self):
        """Execute transformation kernel
        """
        # TODO: should be exposed by the helper
        kernel = self.helper.kernels.get_kernel("transform")
        localshape = self.helper.wg["transform"]
        inshape = self.helper.shape[::-1]
        outshape = self.helper.outshape[::-1]
        globalshape = calc_size(inshape, localshape)
        cval = self.helper.sift.cl_mem["min"].get()[0]
        kernel(
            self._gpu_queue,
            globalshape,
            localshape,
            self._gpu_buffers["input"].data,
            self._gpu_buffers["output"].data,
            self._gpu_buffers["matrix"].data,
            self._gpu_buffers["offset"].data,
            np.int32(inshape[0]),
            np.int32(inshape[1]),
            np.int32(outshape[0]),
            np.int32(outshape[1]),
            cval,
            np.int32(1),
        )  # bilinear interpolation

    def set_reference(self, img, previous=False):
        """Reference for alignment
        """
        if previous:
            img = self._gpu_buffers["input"].data
        self.helper = sift.LinearAlign(img)
        data = np.ascontiguousarray([self.cval], self.siftdtype)
        cpy = pyopencl.enqueue_copy(
            self._gpu_queue, self.helper.sift.cl_mem["min"].data, data
        )
        cpy.wait()

    def get_alignkernel(self):
        """Get transformation from alignment kernel
        """
        return self._transform

    def set_transformkernel(self, transform):
        """Set the transformation kernel according to the alignment kernel
        and adapted transformation
        """
        self._transform.fromtransform(transform)
        self._transform_to_buffer()
