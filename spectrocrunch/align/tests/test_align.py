# -*- coding: utf-8 -*-

import unittest

from .. import alignElastix
from .. import alignSift
from ..alignFFT import alignFFT
from ..alignSimple import alignMin
from ..alignSimple import alignMax
from ..alignSimple import alignCentroid
from ..alignSimple import alignGaussMax
from ..types import transformationType
from .import helper_teststack

import numpy as np


siftctx = alignSift.create_context()


class test_align(unittest.TestCase):

    def compare_relativecof(self, cofs, cofrel, msg=None):
        for i in range(1, cofs.shape[0]):
            cofrelcalc = np.dot(np.linalg.inv(cofs[i-1, ...]), cofs[i, ...])
            np.testing.assert_almost_equal(
                cofrelcalc, cofrel, decimal=1, err_msg=msg)

    def try_alignment(self, alignclass, transfotype, realistic=False, subpixel=True):
        if transfotype == transformationType.translation and\
                alignclass != alignSift.alignSift and\
                alignclass != alignElastix.alignElastix:
            lst = [False, True]
        else:
            lst = [False]

        for vector in lst:
            for transposed in lst:
                if transposed and not vector:
                    continue

                # Prepare dataIO
                inputstack, cofrel, stackdim = helper_teststack.data(
                    transfotype, vector=vector, transposed=transposed, realistic=realistic, subpixel=subpixel)
                outputstack = [np.zeros(1, dtype=np.float32)]*len(inputstack)

                # References
                refdatasetindex = 0
                refimageindex = 0  # len(inputstack)//2

                # Prepare alignment
                o = alignclass(inputstack, None, outputstack, None, None, stackdim=stackdim,
                               overwrite=True, plot=False, transfotype=transfotype)

                # Check alignment
                if vector:
                    if transposed:
                        roi = ((1, -3), (0, 1))
                    else:
                        roi = ((0, 1), (1, -3))
                else:
                    roi = ((1, -3), (1, -3))

                for i in range(4):
                    pad = (i & 1) == 1
                    crop = (i & 2) == 2

                    # Fixed reference
                    msg = "Alignment: Pad = {}, Crop = {}, 1D = {}, transposed = {}, type = {}".format(
                        pad, crop, vector, transposed, "fixed")
                    o.align(refdatasetindex, refimageindex=refimageindex,
                            pad=pad, crop=crop, roi=roi)
                    self.compare_relativecof(o.absolute_cofs(
                        homography=True), cofrel, msg=msg)

                    # Pairwise: align on raw
                    msg = "Alignment: Pad = {}, Crop = {}, 1D = {}, transposed = {}, type = {}".format(
                        pad, crop, vector, transposed, "pairwise/raw")
                    o.align(refdatasetindex, onraw=True,
                            pad=pad, crop=crop, roi=roi)
                    self.compare_relativecof(o.absolute_cofs(
                        homography=True), cofrel, msg=msg)

                    # Pairwise: align on aligned
                    msg = "Alignment: Pad = {}, Crop = {}, 1D = {}, transposed = {}, type = {}".format(
                        pad, crop, vector, transposed, "pairwise")
                    o.align(refdatasetindex, onraw=False,
                            pad=pad, crop=crop, roi=roi)
                    self.compare_relativecof(o.absolute_cofs(
                        homography=True), cofrel, msg=msg)

    @unittest.skip("TODO")
    def test_fft_internals(self):
         # Initialize alignFFT (not important)
        inputstack = [np.zeros((2, 2, 2), dtype=np.float32)]*5
        outputstack = [np.zeros(1, dtype=np.float32)]*5
        o = alignFFT(inputstack, None, outputstack, None, None, stackdim=2,
                     overwrite=True, transfotype=transformationType.similarity)

        # Test Fourier related things
        img = np.abs(np.fft.fft2(np.arange(7*8).reshape((7, 8))))
        img2 = np.fft.ifftshift(np.fft.fftshift(img))
        np.testing.assert_allclose(img, img2)

        # Test log-polar
        mx = 10
        my = 10
        nx = 501
        ny = 401
        dx = 2*mx/(nx-1.)
        dy = 2*my/(ny-1.)
        xv, yv = np.meshgrid(np.linspace(-mx, mx, nx),
                             np.linspace(-my, my, ny))

        sx = 2.
        sy = 1.
        angle = 0.
        data = [(0., 0., sx, sy, angle, 1000.)]
        fixed = helper_teststack.gettransformedimage(
            xv, yv, data, angle=True).reshape(xv.shape)

        angle = -2*np.pi/180
        scale = 1.3  # 0.9
        sx *= scale
        sy *= scale
        data = [(0., 0., sx, sy, angle, 1000.)]
        moving = helper_teststack.gettransformedimage(
            xv, yv, data, angle=True).reshape(xv.shape)

        a = scale*np.cos(angle)
        b = scale*np.sin(angle)
        R = np.array([[a, -b, 0], [b, a, 0], [0, 0, 1]])

        T = np.array([[1, 0, mx/dx], [0, 1, my/dy], [0, 0, 1]])
        Tinv = np.array([[1, 0, -mx/dx], [0, 1, -my/dy], [0, 0, 1]])

        M = np.dot(T, np.dot(R, Tinv))

        o.set_reference(fixed)
        aligned = o.execute_alignkernel(moving)

        np.testing.assert_almost_equal(
            M, o._transform.getnumpyhomography(), decimal=1)

    @unittest.skipIf(alignElastix.sitk is None, "SimpleElastix is not installed")
    def test_elastix(self):
        types = [transformationType.translation]
        for t in types:
            self.try_alignment(alignElastix.alignElastix, t)

    @unittest.skipIf(alignSift.pyopencl is None, "pyopencl is not installed")
    @unittest.skipIf(siftctx is None, "no pyopencl device available")
    @unittest.skipIf(True, "temporary disable")
    def test_sift(self):
        types = [transformationType.translation, transformationType.rigid,
                 transformationType.similarity, transformationType.affine]
        # TODO: the others are not that precise
        types = [transformationType.translation]
        for t in types:
            self.try_alignment(alignSift.alignSift, t)

    def test_fft(self):
        types = [transformationType.translation]
        for t in types:
            self.try_alignment(alignFFT, t)

    def test_min(self):
        self.try_alignment(
            alignMin, transformationType.translation, realistic=True, subpixel=False)

    def test_max(self):
        self.try_alignment(
            alignMax, transformationType.translation, subpixel=False)

    def test_centroid(self):
        self.try_alignment(
            alignCentroid, transformationType.translation, subpixel=False)

    def test_gaussmax(self):
        self.try_alignment(
            alignGaussMax, transformationType.translation, subpixel=False)


def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_align("test_fft_internals"))
    testSuite.addTest(test_align("test_min"))
    testSuite.addTest(test_align("test_max"))
    testSuite.addTest(test_align("test_centroid"))
    testSuite.addTest(test_align("test_gaussmax"))
    testSuite.addTest(test_align("test_fft"))
    testSuite.addTest(test_align("test_elastix"))
    testSuite.addTest(test_align("test_sift"))
    return testSuite


if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
