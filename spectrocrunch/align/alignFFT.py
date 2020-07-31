# -*- coding: utf-8 -*-

import numpy as np
import scipy.ndimage.interpolation

from .align import align
from .types import transformationType
from ..math import center
from ..math import ft


class alignFFT(align):
    def __init__(self, *args, **kwargs):
        super(alignFFT, self).__init__(*args, **kwargs)

        # FT's of images
        self.fixedft = None
        self.movingft = None

        # change of reference frame
        self._transform = self.defaulttransform()

        # subpixel sampling factor > 1
        self.sampling = 20  # shift +/- 0.05
        # self.sampling = 100 # shift +/- 0.01

        # maximum phase correlation: max, fit, centroid
        self.maxmethod = "centroid"

    def execute_transformkernel(self, img):
        """Transform image according with the transformation kernel
        """
        return self._transform.transformimage(img)

    def ifft_interpolate(self, imgft, ROIoffset, ROIsize):
        """Sub-region inverse Fourier transform with subpixel interpolation
             Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup,
             "Efficient subpixel image registration algorithms,"
             Optics Letters 33, 156-158 (2008).

        Args:
            imgft: FFT of images
            ROIoffset: sub-pixel grid (super-space)
            ROIsize: sub-pixel grid (super-space)

        Returns:
            np.array: inverse FFT on sub-pixel grid
        """

        if imgft.size in imgft.shape:
            # Super space
            x0 = -ROIoffset
            x1 = x0 + ROIsize - 1

            # Frequencies corresponding to super space
            nu = imgft.size
            u = ft.fftfreq(nu, d=self.sampling)

            return ft.ifft(imgft, x0=x0, x1=x1, u=u)
        else:
            # Super space
            x0 = -ROIoffset[1]
            x1 = x0 + ROIsize[1] - 1
            y0 = -ROIoffset[0]
            y1 = y0 + ROIsize[0] - 1

            # Frequencies corresponding to super space
            nu = imgft.shape[1]
            nv = imgft.shape[0]
            u = ft.fftfreq(nu, d=self.sampling)
            v = ft.fftfreq(nv, d=self.sampling)

            return ft.ifft2(imgft, x0=x0, x1=x1, u=u, y0=y0, y1=y1, v=v)

    def logpolar(self, imgft, scale_estimate=None):
        """Log-polar coordinate transformation of the modulus of the given Fourier transform
        """

        # Filtered modulus of the ft
        imgftmod = self.ftmodulus(imgft)
        # self.plot2(ft.fftshift(imgftmod),fourierlines=True,title="|FFT|")
        imgftmod = self.highpassfilter(imgftmod)
        # self.plot2(ft.fftshift(imgftmod),fourierlines=True,title="|FFT| with HP filter")
        imgftmod = ft.fftshift(imgftmod)

        # Size and frequency limits
        s = np.array(imgft.shape)
        sy, sx = imgft.shape
        pmin, pmax = ft.fft_freqind(s)

        # Log-polar grid for interpolation
        if scale_estimate:
            sclmin = scale_estimate

        rmax = np.maximum(s - 1 - pmax, pmax)
        rmax = np.linalg.norm(rmax)

        dt = np.pi / (sy - 1)
        b = 10 ** (np.log10(rmax) / (sx - 1))
        r = b ** np.arange(sx)
        t = dt * np.arange(sy)

        # Log-polar to Cartesian
        x = np.outer(np.cos(t), r) - pmin[1]
        y = np.outer(np.sin(t), r) - pmin[0]

        # Interpolation
        out = np.empty_like(x)
        scipy.ndimage.interpolation.map_coordinates(imgftmod, [y, x], output=out)

        return out, b, dt

    def plot2(self, img, index=2, fourierlines=False, title="none"):
        """For testing purposes
        """
        import matplotlib.pyplot as plt

        fig = plt.figure(index)
        if 0:
            import mpl_toolkits.mplot3d as p3

            ax = p3.Axes3D(fig)
            yy, xx = np.indices(img.shape)
            ax.plot_surface(yy, xx, img)
        else:
            ax = fig.add_subplot(111)

            if img.size in img.shape:
                x = ft.fftfreq(len(img), 1, centered=True)
                ax.plot(x, img)
            else:
                ax.imshow(img, origin="lower", interpolation="nearest")
                if fourierlines:
                    s = np.asarray(img.shape)
                    fmin, fmax = ft.fft_freqind(s)
                    cen = -fmin
                    ax.plot([cen[1], cen[1]], [0, s[0] - 1])
                    ax.plot([0, s[1] - 1], [cen[0], cen[0]])

        plt.title(title)
        plt.pause(0.1)
        raw_input("press enter...")

    def movingFTdummy(self):
        """FT of moving image with a known shift w.r.t. the reference
        """
        d0, d1 = self.fixedft.shape
        v0, v1 = (3.48574, 8.73837)
        f0 = ft.ifftshift(np.arange(-d0 // 2, d0 // 2))
        f1 = ft.ifftshift(np.arange(-d1 // 2, d1 // 2))
        m1, m0 = np.meshgrid(f1, f0)
        e0 = np.exp(-2j * np.pi * v0 * m0 / float(d0))
        e1 = np.exp(-2j * np.pi * v1 * m1 / float(d1))
        e = e0 * e1
        self.movingft = self.fixedft * e

    def fft_handle_missing(self, img):
        """FFT with handling of missing data
        """
        if self.cval != 0:
            if self.cval is np.nan:
                missing = np.isnan(img)
            else:
                missing = img == self.cval
            bmissing = np.any(missing)
        else:
            bmissing = False

        if bmissing:
            img2 = img.copy()
            img2[missing] = 0
        else:
            img2 = img

        if img2.size in img.shape:
            return ft.fft(img2.flatten())
        else:
            return ft.fft2(img2)

    def ftmodulus(self, imgft):
        imgabs = np.abs(imgft)
        idx = imgabs < 1.0e-20
        if idx.any():
            imgabs[idx] = 1.0
        return imgabs

    def highpassfilter(self, imgft):
        """High-pass filter in Fourier domain
        """
        nu = imgft.shape[1]
        nv = imgft.shape[0]
        h = np.outer(
            np.cos(ft.fftfreq(nv, d=1 / np.pi)), np.cos(ft.fftfreq(nu, d=1 / np.pi))
        )
        return imgft * (1 - h) / (2 - h)

    def cc_maximum(self, cross_correlation, maxmethod="max"):
        data = self.ftmodulus(cross_correlation)
        # data = cross_correlation.real

        if maxmethod == "centroid":
            shift = center.fcentroid(data)
        elif maxmethod == "fit":
            shift = center.fgaussmax(data)
        else:
            shift = center.fmax(data)

        return shift

    def determine_shift(self, img1ft, img2ft):
        """Determine shift between images using their Fourier transforms
        """
        is1D = img1ft.size in img1ft.shape

        # Calculate shift without subpixel precision
        # Ideally this should be a delta function: real(F.G*/(|F|.|G|))
        # In reality this doesn't seem to work, somehow |F.G*| is better, no idea why
        image_product = img1ft * img2ft.conj()
        # image_product /= self.ftmodulus(image_product)
        if is1D:
            cross_correlation = ft.ifft(image_product)
        else:
            cross_correlation = ft.ifft2(image_product)
        shift = self.cc_maximum(cross_correlation)
        # self.plot2(ft.fftshift(self.ftmodulus(cross_correlation)),fourierlines=True,title='|cross correlation|')

        # Shift indices = [0,...,imax,imin,...,-1]
        if is1D:
            _, imax = ft.fft_freqind(cross_correlation.size)
            if shift > imax:
                shift -= cross_correlation.size
        else:
            s = np.array(cross_correlation.shape)
            _, imax = ft.fft_freqind(s)
            shift[shift > imax] -= s[shift > imax]

        # print("Initial shift {}".format(shift))

        # Calculate shift with subpixel precision by interpolating
        # the cross-correlation around the maximum (or center-of-mass).
        if self.sampling > 1:
            # real-space: n, d=1, shift=s
            # super-space: n.sampling, d=1, shift=s.sampling
            if is1D:
                ROIsize = 4 * self.sampling
            else:
                # ROI in super-space = 3x3 pixels in real-space
                ROIsize = self.sampling * np.array((3, 3))
            _, ksampled = ft.fft_freqind(ROIsize)  # ROI center

            # ROI = [0,1,...,ROIsize-1] - ksampled + shift*self.sampling  (so that maximum in the middle)
            ROIoffset = self.dtype(ksampled) - shift * self.sampling
            cross_correlation = self.ifft_interpolate(image_product, ROIoffset, ROIsize)
            # cross_correlation /= img2ft.shape[0]*img2ft.shape[1] * self.sampling ** 2
            # self.plot2(self.ftmodulus(cross_correlation),title='|cross correlation super space|')

            # Get fit, maximum, centroid, ...
            shiftsampled = self.cc_maximum(cross_correlation, maxmethod=self.maxmethod)

            # Shift from super space to real space
            shift = (shiftsampled - ROIoffset) / self.dtype(self.sampling)

        # import skimage.feature
        # shft,tmp1,tmp2 = skimage.feature.register_translation(img1ft, img2ft, upsample_factor=self.sampling,space="fourier")
        # print(shft,shift)
        # assert(np.array_equal(shft,shift))

        # print("Final shift {}".format(shift))

        return shift

    def execute_alignkernel(self, img):
        """Align image on reference"""

        self.movingft = self.fft_handle_missing(img)

        if self.movingft.size in self.movingft.shape:
            if self.transfotype != transformationType.translation:
                raise ValueError("Only translations apply to 1D data.")
            else:
                shift = self.determine_shift(self.movingft, self.fixedft)
                if img.shape[0] == 1:
                    self._transform.settranslation(shift, 0)
                else:
                    self._transform.settranslation(0, shift)
        else:
            if self.transfotype == transformationType.translation:
                shift = self.determine_shift(self.movingft, self.fixedft)
                self._transform.settranslation(shift[1], shift[0])
            elif (
                self.transfotype == transformationType.similarity
                or self.transfotype == transformationType.rigid
            ):
                # raise NotImplementedError("alignFFT doesn't support this type transformation.")

                # Fourier transforms in log-polar coordinates
                fixedft_lp, _, _ = self.logpolar(self.fixedft)
                movingft_lp, b, dt = self.logpolar(self.movingft)
                self.plot2(fixedft_lp, title="Fixed LP")

                # Determine shift between the log-polar images
                fixedft_lp = ft.fft2(fixedft_lp)
                movingft_lp = ft.fft2(movingft_lp)
                shift = self.determine_shift(movingft_lp, fixedft_lp)

                # Convert shift from pixels to scale and angle
                if self.transfotype == transformationType.rigid:
                    scale = 1
                else:
                    scale = b ** (-shift[1])
                angle = -shift[0] * dt

                print(b ** (np.arange(-shift[1] - 2, -shift[1] + 5)))
                print(-np.arange(shift[0] - 2, shift[0] + 5) * dt * 180 / np.pi)
                print(scale, angle * 180 / np.pi)
                # scale = 1.3
                # angle = 2*np.pi/180
                # print(scale,angle*180/np.pi)

                # Linear part of the similarity transform
                a = scale * np.cos(angle)
                b = scale * np.sin(angle)
                self._transform.setaffine(np.array([[a, b, 0], [-b, a, 0]]))

                # Remove rotation and scaling from moving image
                moving = self.execute_transformkernel(img)
                movingft = self.fft_handle_missing(moving)

                # Determine shift (this doesn't work well, maybe because of missing data?)
                shift = self.determine_shift(movingft, self.fixedft)
                shift = [-48.0, -84.0]
                self._transform.settranslation(shift[1], shift[0])

                moving = self.execute_transformkernel(img)
                tmp = np.zeros(self.fixedft.shape + (3,), dtype=moving.dtype)
                tmp[..., 0] = ft.ifft2(self.fixedft).real
                tmp[..., 1] = moving
                tmp /= np.nanmax(tmp)
                self.plot2(tmp, title="RGB")

            else:
                raise NotImplementedError(
                    "alignFFT doesn't support this type transformation."
                )

        return self.execute_transformkernel(img)

    def set_reference(self, img, previous=False):
        """Reference for alignment
        """
        if previous:
            self.fixedft = self.movingft
        else:
            self.fixedft = self.fft_handle_missing(img)

    def get_alignkernel(self):
        """Get transformation
        """
        return self._transform

    def set_transformkernel(self, transfo):
        """Set transformation
        """
        self._transform.fromtransform(transfo)
