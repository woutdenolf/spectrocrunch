# -*- coding: utf-8 -*-
#
#   Copyright (C) 2016 European Synchrotron Radiation Facility, Grenoble, France
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

from ..utils.Enum import Enum
fftConvention = Enum(['numpy','idl'])

FFT_FREQ_CONVENTION = fftConvention.numpy
# n: number of data points in real space
# d: pixel spacing in real space
#    Convention numpy:
#     frequency = [0,...,k,-k-1,...-1]/(n.d)  (n is even, k = n/2-1  , -k-1 = -n//2, nfreq = k+1 + k+1 = n )
#                 [0,...,k,-k,...,-1]/(n.d)   (n is odd,  k = (n-1)/2, -k   = -n//2, nfreq = k+1 + k = n )
#                 k = (n+1)//2-1
#            
#    Convention IDL:
#     frequency = [0,...,k,-k+1,...-1]/(n.d)  (n is even, k = n/2    , -k+1 = -n//2+1, nfreq = k+1 + k-1 = n )
#                 [0,...,k,-k,...,-1]/(n.d)   (n is odd,  k = (n-1)/2, -k   = -n//2  , nfreq = k+1 + k = n )
#                 k = (n+1)//2-(n mod 2)

FFT_NORM_CONVENTION = fftConvention.numpy
# n: number of data points in real space
# Convention numpy:
#   FT  = sum(exp(-2.pi.i.u.x))
#   IFT = sum(exp(2.pi.i.u.x))/n
# Convention IDL ():
#   FT  = sum(exp(-2.pi.i.u.x))/n
#   IFT = sum(exp(2.pi.i.u.x))

import numpy as np

def fft_freqind(n,freqconvention=FFT_FREQ_CONVENTION):
    """Calculate frequency range of an fft (multiplied by n.d)

        Args:
            n (int|np.array): number of data points in real space
            freqconvention (Optional(bool)): even frequency convention 1=numpy, 2=IDL
    """
    # frequencies = [0,...,imax,imin,...,-1]/(n.d)
    if freqconvention==fftConvention.idl:
        imin = -(n//2)+(1-(n%2))
        imax = (n+1)//2-(n%2)
    else:
        imin = -(n//2)
        imax = (n+1)//2-1
    return imin,imax

def fftfreq(n,d=1,centered=False,freqconvention=FFT_FREQ_CONVENTION):
    """Fourier space with zero frequency first

        Args:
            n (int): number of real space data points
            d (num): real space data point spacing
            freqconvention (Optional(bool)): even frequency convention 1=numpy, 2=IDL
    """
    imin,imax = fft_freqind(n,freqconvention=freqconvention)
    if centered:
        freq = np.arange(imin, imax+1, dtype=int)
    else:
        freq = np.empty(n, dtype=int)
        npos = imax+1
        freq[:npos] = np.arange(npos, dtype=int)
        freq[npos:] = np.arange(imin, 0, dtype=int)
    return freq/float(n*d)

def fftshift(sigft,freqconvention=FFT_FREQ_CONVENTION):
    """Shift zero frequency to the middle

        Args:
            sigft (np.array): signal in Fourier space
            centered (Optional(bool)): zero frequency in the middle
            freqconvention (Optional(bool)): even frequency convention 1=numpy, 2=IDL
    """
    dim = np.array(sigft.shape)
    _, imax = fft_freqind(dim,freqconvention=freqconvention)
    npos = imax+1
    out = sigft.copy()
    for k in range(len(dim)):
        ind = np.empty(dim[k], int)
        off = dim[k]-npos[k]
        ind[:off] = np.arange(npos[k], dim[k], dtype=int)
        ind[off:] = np.arange(npos[k], dtype=int)
        out = np.take(out, ind, axis=k)
    return out

def ifftshift(sigft,freqconvention=FFT_FREQ_CONVENTION):
    """Shift zero frequency to zero

        Args:
            sigft (np.array): signal in Fourier space
            freqconvention (Optional(bool)): even frequency convention 1=numpy, 2=IDL
    """
    dim = np.array(sigft.shape)
    imin, _ = fft_freqind(dim,freqconvention=freqconvention)
    npos = -imin
    out = sigft.copy()
    for k in range(len(dim)):
        ind = np.empty(dim[k], int)
        off = dim[k]-npos[k]
        ind[:off] = np.arange(npos[k], dim[k], dtype=int)
        ind[off:] = np.arange(npos[k], dtype=int)
        out = np.take(out, ind, axis=k)
    return out

def _realspace(N,dx,x0,x1):
    """Real space data points must be equally spaced (but they might be a subset)

        Args:
            N (int): number of real or Fourier space data points
            dx (num): real space increment
            x0 (int): real space start integer index
            x1 (int): real space end integer index
    """
    if x1 is None:
        x1 = N-1
    if x0 > x1:
        raise ValueError("Wrong real space coordinates ({},{})".format(x0,x1))
    return np.arange(x0,x1+1)*dx

def _dft(f,dx=1,x0=0,x1=None,u=[],centered=False,inverse=False,normconvention=FFT_NORM_CONVENTION):
    """Fourier transform with fixed frequencies

        Args:
            f (np.array): function in real or Fourier space
            dx (Optional(num)): real space data point spacing
            x0 (Optional(num)): real space start index
            x1 (Optional(num)): real space end index
            u (Optional(np.array)): Fourier space
            centered (Optional(bool)): zero frequency in the middle
            inverse (Optional(bool)): inverse Fourier transform
            normconvention (Optional(bool)): fft normalization 1=numpy, 2=IDL
    """

    if dx==1 and x0==0 and x1 is None and len(u)==0 and not centered:
        if inverse:
            ret = np.fft.ifft(f)
            if normconvention==fftConvention.idl:
                ret *= len(f)
        else:
            ret = np.fft.fft(f)
            if normconvention==fftConvention.idl:
                ret /= len(f)
    else:
        # Real space
        x = _realspace(len(f),dx,x0,x1)

        # Fourier space
        if len(u)==0:
            u = fftfreq(len(f),dx,centered=centered)

        # Check dimensions
        if inverse:
            if len(u)!=len(f):
                raise ValueError("Number of frequencies should be equal to the number of data points in Fourier space")
            c = 2j * np.pi
            ret = np.exp(c*np.outer(x,u)).dot(f) # nx x nu x nu
        else:
            if len(x)!=len(f):
                raise ValueError("Number of times should be equal to the number of data points in real space")
            c = -2j * np.pi
            ret = np.exp(c*np.outer(u,x)).dot(f) # nu x nx x nx

        # Normalization:
        if (normconvention==fftConvention.idl) ^ inverse:
            ret /= len(u)

    return ret

def _dft2(f,dx=1,x0=0,x1=None,u=[],\
            dy=1,y0=0,y1=None,v=[],\
            centered=False,inverse=False,\
            normconvention=FFT_NORM_CONVENTION):
    """Fourier transform with fixed frequencies
        Sub-region inverse Fourier transform with subpixel interpolation using the matrix for of the 2D-DFT
            Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup,
            "Efficient subpixel image registration algorithms,"
            Optics Letters 33, 156-158 (2008).

        Args:
            f (np.array): function in real or Fourier space
            dx (Optional(num)): real space data point spacing
            x0 (Optional(num)): real space start index
            x1 (Optional(num)): real space end index
            u (Optional(np.array)): Fourier space
            dy (Optional(num)): real space data point spacing
            y0 (Optional(num)): real space start index
            y1 (Optional(num)): real space end index
            v (Optional(np.array)): Fourier space
            centered (Optional(bool)): zero frequency in the middle
            inverse (Optional(bool)): inverse Fourier transform
            normconvention (Optional(bool)): fft normalization 1=numpy, 2=IDL
    """

    if dx==1 and x0==0 and x1 is None and len(u)==0 and\
       dy==1 and y0==0 and y1 is None and len(v)==0 and not centered:
        if inverse:
            ret = np.fft.ifft2(f)
            if normconvention==fftConvention.idl:
                ret *= f.shape[0]*f.shape[1]
        else:
            ret = np.fft.fft2(f)
            if normconvention==fftConvention.idl:
                ret /= f.shape[0]*f.shape[1]
    else:
        # Real space
        x = _realspace(f.shape[1],dx,x0,x1)
        y = _realspace(f.shape[0],dy,y0,y1)

        # Fourier space
        if len(u)==0:
            u = fftfreq(f.shape[1],dx,centered=centered)
        if len(v)==0:
            v = fftfreq(f.shape[0],dy,centered=centered)

        # DFT (forward or backward)
        if inverse:
            if len(u)!=f.shape[1] or len(v)!=f.shape[0]:
                raise ValueError("Number of frequencies should be equal to the number of data points in Fourier space")
            c = 2j * np.pi
            col_kernel = np.exp(c * u[:, None].dot(x[None, :])) # nu x nx
            row_kernel = np.exp(c * y[:, None].dot(v[None, :])) # ny x nv
            # ny x nv . nv x nu . nu x nx
        else:
            if len(x)!=f.shape[1] or len(y)!=f.shape[0]:
                raise ValueError("Number of times should be equal to the number of data points in Fourier space")
            c = -2j * np.pi
            col_kernel = np.exp(c * x[:, None].dot(u[None, :])) # nx x nu
            row_kernel = np.exp(c * v[:, None].dot(y[None, :])) # nv x ny
            # nv x ny . ny x nx . nx x nu
        ret = (row_kernel.dot(f).dot(col_kernel)) 

        # Normalization:
        if (normconvention==fftConvention.idl) ^ inverse:
            ret /= (len(u)*len(v))

    return ret

def fft(f,**kwargs):
    """Fourier transform with default frequencies

        Args:
            f (np.array): function in real space
    """
    return _dft(f,inverse=False,**kwargs)

def ifft(f,**kwargs):
    """Inverse Fourier transform with default frequencies

        Args:
            f (np.array): function in Fourier space
    """
    return _dft(f,inverse=True,**kwargs)

def fft2(f,**kwargs):
    """Fourier transform with default frequencies

        Args:
            f (np.array): function in real space
    """
    return _dft2(f,inverse=False,**kwargs)

def ifft2(f,**kwargs):
    """Inverse Fourier transform with default frequencies

        Args:
            f (np.array): function in real space
    """
    return _dft2(f,inverse=True,**kwargs)


