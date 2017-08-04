# -*- coding: utf-8 -*-
#
#   Copyright (C) 2017 European Synchrotron Radiation Facility, Grenoble, France
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

from .. import xiaedf

import numpy as np

def random(a,b,n):
    return a+(b-a)*np.random.random(n)

def peak(x, p):
    h, mu, sig = p
    return h*np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def ctsround(x,stattype):
    return stattype(round(x))

def data(nspec,nchan,ndet):
    datatype = np.float
    stattype = np.int

    # Detector characteristics
    solidangle = np.linspace(1,1.5,ndet,dtype=datatype)
    DT = 2*np.arange(ndet)
        
    # Peak characteristics
    npeaks = 10
    x = np.arange(nchan,dtype=datatype)
    p = zip(random(1000,2000,npeaks),\
            random(0,nchan,npeaks),\
            random(20,30,npeaks))

    # Generate data
    dataorg = np.zeros((nspec,nchan,ndet),dtype=datatype)
    data = np.zeros((nspec,nchan,ndet),dtype=datatype)
    stats = np.zeros((nspec,xiaedf.xiadata.NSTATS,ndet),dtype=stattype)
    for i in range(nspec):
        spectrum = np.zeros(nchan)
        for k in range(npeaks):
            spectrum += peak(x,p[k])

        for j in range(ndet):
            rawspectrum = spectrum*solidangle[j]
            rawspectrum = rawspectrum.astype(stattype)
            rawspectrum += np.random.poisson(rawspectrum)
            ICR = rawspectrum.sum()

            OCR = ICR*(1-j/50.)# DT = 2*j %
            spectrumj = rawspectrum*OCR/ICR
            spectrumj = spectrumj.astype(stattype)
            OCR = spectrumj.sum()

            data[i,:,j] = spectrumj

            dataorg[i,:,j] = spectrumj*(datatype(ICR)/OCR)

            stats[i,xiaedf.xiadata.STDET,j] = j
            stats[i,xiaedf.xiadata.STEVT,j] = ICR # % Not sure
            stats[i,xiaedf.xiadata.STICR,j] = ICR
            stats[i,xiaedf.xiadata.STOCR,j] = OCR
            stats[i,xiaedf.xiadata.STDT,j] = ctsround(100-OCR*100./ICR,stattype) # %
            stats[i,xiaedf.xiadata.STLT,j] = ctsround(OCR*1000./ICR,stattype) # 1000 msec RT

    return dataorg,data,stats

