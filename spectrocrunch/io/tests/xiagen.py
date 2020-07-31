# -*- coding: utf-8 -*-

from .. import xiaedf

import numpy as np


def random(a, b, n):
    return a + (b - a) * np.random.random(n)


def peak(x, p):
    h, mu, sig = p
    return h * np.exp(-np.power(x - mu, 2.0) / (2 * np.power(sig, 2.0)))


def ctsround(x, stattype):
    return stattype(round(x))


def data(nspec, nchan, ndet, concentration=None, flux=None):
    datatype = np.int32
    stattype = np.int32
    cortype = np.float64

    # Detector characteristics
    solidangle = np.linspace(1, 1.5, ndet, dtype=cortype)
    DT = 2 * np.arange(ndet)

    # Peak characteristics
    npeaks = 10
    x = np.arange(nchan, dtype=datatype)
    p = zip(
        random(1000, 2000, npeaks), random(0, nchan, npeaks), random(20, 30, npeaks)
    )

    # Spectrum
    spectrum0 = np.zeros(nchan, dtype=cortype)
    for k in p:
        spectrum0 += peak(x, k)

    # Pixel intensity
    if concentration is None:
        concentration = np.ones(nspec, dtype=cortype)

    # Pixel flux
    if flux is None:
        flux = np.ones(nspec, dtype=cortype)

    # Generate data
    dataorg = np.zeros((nspec, nchan, ndet, 3), dtype=cortype)
    data = np.zeros((nspec, nchan, ndet), dtype=datatype)
    stats = np.zeros((nspec, xiaedf.xiadata.NSTATS, ndet), dtype=stattype)
    for i in range(nspec):
        spectrum = spectrum0 * concentration[i] * flux[i]  # cortype

        for j in range(ndet):
            rawspectrum = spectrum * solidangle[j]
            rawspectrum = np.random.poisson(rawspectrum)
            ICR = rawspectrum.sum()

            OCR = ICR * (1 - j / 50.0)  # DT = 2*j %
            spectrumj = rawspectrum * OCR / ICR
            OCR = spectrumj.sum()

            spectrumj = spectrumj.astype(datatype)
            ICR = ICR.astype(stattype)
            OCR = OCR.astype(stattype)

            data[i, :, j] = spectrumj

            dataorg[i, :, j, 0] = spectrumj * (cortype(ICR) / cortype(OCR))
            dataorg[i, :, j, 1] = spectrumj / cortype(flux[i])
            dataorg[i, :, j, 2] = (
                spectrumj * ((cortype(ICR) / cortype(OCR))) / cortype(flux[i])
            )

            stats[i, xiaedf.xiadata.STDET, j] = j
            stats[i, xiaedf.xiadata.STEVT, j] = ICR  # % Not sure
            stats[i, xiaedf.xiadata.STICR, j] = ICR
            stats[i, xiaedf.xiadata.STOCR, j] = OCR
            stats[i, xiaedf.xiadata.STDT, j] = ctsround(
                100 - OCR * 100.0 / ICR, stattype
            )  # %
            stats[i, xiaedf.xiadata.STLT, j] = ctsround(
                OCR * 1000.0 / ICR, stattype
            )  # 1000 msec RT

    return dataorg, data, stats
