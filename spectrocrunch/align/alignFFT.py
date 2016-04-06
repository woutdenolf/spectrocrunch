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
import numpy as np
import scipy.ndimage

class alignFFT(align):

    def __init__(self,*args,**kwargs):
        super(alignFFT,self).__init__(*args,**kwargs)

        # FT's of images
        self.fixedft = None
        self.movingft = None
        
        # change of reference frame
        self.offset = self.idoffset.copy()
        self.linear = self.idlinear.copy()

        # subpixel sampling factor > 1
        self.sampling = 20

    def execute_transformkernel(self,img):
        """Transform image according with the transformation kernel
        """
        return self.execute_transform_nokernel(img,self.offset,self.linear)
 
    def ifft_interpolate(self,imgft,ROIoffset,ROIsize):
        """Sub-region inverse Fourier transform with subpixel interpolation using the matrix for of the 2D-DFT
        """
        nu = imgft.shape[1]
        nv = imgft.shape[0]
        u = np.fft.fftfreq(nu,self.sampling)[:, None] # nu x 1
        v = np.fft.fftfreq(nv,self.sampling)[None, :] # 1 x nv

        x = np.arange(ROIsize[1])[None, :] - ROIoffset[1] # 1 x nxsub
        y = np.arange(ROIsize[0])[:, None] - ROIoffset[0] # nysub x 1
        
        c = 2j * np.pi
        col_kernel = np.exp(c * u.dot(x)) # nu x nxsub
        row_kernel = np.exp(c * y.dot(v)) # nysub x nv

        # nysub x nv . nv x nu . nu x nxsub
        return (row_kernel.dot(imgft).dot(col_kernel))

    def plot2(self,img):
        import pylab
        fig = pylab.figure(2)
        if 0:
            import mpl_toolkits.mplot3d as p3
            ax = p3.Axes3D(fig)
            xx, yy = np.mgrid[0:img.shape[0],0:img.shape[1]]
            ax.plot_surface(xx,yy,img)
        else:
            pylab.imshow(img,origin='lower',interpolation='nearest')
        pylab.pause(1)

    def movingFTdummy(self):
        """FT of moving image with a known shift w.r.t. the reference
        """
        d0, d1 = self.fixedft.shape
        v0, v1 = (3.48574,8.73837)
        f0 = np.fft.ifftshift(np.arange(-d0 // 2, d0 // 2))
        f1 = np.fft.ifftshift(np.arange(-d1 // 2, d1 // 2))
        m1, m0 = np.meshgrid(f1, f0)
        e0 = np.exp(-2j * np.pi * v0 * m0 / float(d0))
        e1 = np.exp(-2j * np.pi * v1 * m1 / float(d1))
        e = e0 * e1
        self.movingft = self.fixedft * e

    def fft_handle_missing(self,img):
        """FFT with handling of missing data
        """
        if self.cval != 0:
            if self.cval is np.nan:
                missing = np.isnan(img)
            else:
                missing = img==self.cval
            bmissing = np.any(missing)
        else:
            bmissing = False

        if bmissing:
            img2 = img.copy()
            img2[missing] = 0
        else:
            img2 = img

        return np.fft.fft2(img2)

    def execute_alignkernel(self,img):
        """Align image on reference

           [1]  Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup,
                "Efficient subpixel image registration algorithms,"
                Optics Letters 33, 156-158 (2008).
        """
        self.movingft = self.fft_handle_missing(img)
        #self.movingFTdummy()
        shape = np.array(self.movingft.shape)

        # Moduli (not used anymore ...)
        #absfixed = abs(self.fixedft)
        #absmoving = abs(self.movingft)
        #idx = absfixed < 1.0e-20
        #if idx.any():
        #    absfixed[idx] = 1.0
        #idx = absmoving < 1.0e-20
        #if idx.any():
        #    absmoving[idx] = 1.0

        # Calculate shift without subpixel precision
        # n: size in real space
        # d: pixel spacing in real space
        #
        # Convention 1 (np.fft.fftfreq, np.fft.fftshift):
        #   frequency = [0,...,k,-k-1,...-1]/(n.d)  (n is even, k = n/2-1  , nfreq = k+1 + k+1 = n )
        #               [0,...,k,-k,...,-1]/(n.d)   (n is odd,  k = (n-1)/2, nfreq = k+1 + k = n )
        #               k = (n+1)//2-1
        #
        # Convention 2 (IDL):
        #   frequency = [0,...,k,-k+1,...-1]/(n.d)  (n is even, k = n/2    , nfreq = k+1 + k-1 = n )
        #               [0,...,k,-k,...,-1]/(n.d)   (n is odd,  k = (n-1)/2, nfreq = k+1 + k = n )
        #               k = (n+1)//2-(n mod 2)
        #
        # We will adopt convention 1. Note that the maximal shift that can be registered is k.
        #
        image_product = self.fixedft * self.movingft.conj() # no need to normalize by abs(self.fixedft).abs(self.movingft) ?

        cross_correlation = np.fft.ifft2(image_product)
        shift = np.array(np.unravel_index(np.argmax(np.abs(cross_correlation)),shape))
        k = (shape+1)//2-1
        shift[shift > k] -= shape[shift > k]

        # Calculate shift with subpixel precision by interpolating
        # the cross-correlation around the maximum (or center-of-mass).
        #
        if self.sampling>1:
            ROIsize = np.ceil(self.sampling * np.array((1.5,1.5))).astype(np.int)
            ksampled = (ROIsize+1)//2-1 
            ROIoffset = ksampled - shift*self.sampling
            cross_correlation = self.ifft_interpolate(image_product,ROIoffset,ROIsize)
            #cross_correlation /= shape[0]*shape[1] * self.sampling ** 2
            #shiftsampled = np.array(np.unravel_index(np.argmax(np.abs(cross_correlation)),ROIsize))
            shiftsampled = np.array(scipy.ndimage.measurements.center_of_mass(np.abs(cross_correlation)),dtype=self.dtype)

            shiftsampled -= ksampled
            shift = shift + shiftsampled / self.dtype(self.sampling)
            #self.plot2(np.abs(cross_correlation))

        #shft,tmp1,tmp2 = skimage.feature.register_translation(self.fixedft, self.movingft, upsample_factor=self.sampling,space="fourier")

        self.offset[:] = -shift
        return self.execute_transformkernel(img)

    def set_reference(self,img,previous=False):
        """Reference for alignment
        """
        if previous:
            self.fixedft = self.movingft
        else:
            self.fixedft = self.fft_handle_missing(img)

    def get_transformation(self):
        """Get transformation
        """
        return self.offset

    def set_transformation(self,offset,changed):
        """Set transformation
        """
        if changed:
            self.offset[:] = offset


