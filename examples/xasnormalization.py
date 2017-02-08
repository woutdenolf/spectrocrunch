# XRF
#  f: identity
#  f': 1
#  f^-1: 1
#  mu(fluo).rho.d = f(I/I0)
#  I = f^-1(m).I0   (m = mu(fluo).rho.d)
# 
# Transmission:
#  f: -ln(x)
#  f': -1/x
#  f^-1: exp(-x)
#  mu(mass).rho.d = f(I/I0)
#  I = f^-1(m).I0   (m = mu(mass).rho.d)
# 
# d(u(x)/v(x)) = u'(x)/v(x) - v'(x)u(x)/[v(x)]^2
# 
# 
# XAS1:
#     XAS1 = f(sum(I)/sum(I0))
#     VARXAS = sumj[ (dXAS/dIj)^2.VARIj + (dXAS/dI0j)^2.VARI0j ]
#
#     dXAS1/dIj  =  f' / sum(I0)
#     dXAS1/dI0j = -f' . sum(I)/[sum(I0)]^2
#
#     Transmission:
#         f' = -sum(I0) / sum(I)
#         dXAS1/dIj  = -1 / sum(I)
#         dXAS1/dI0j = 1 / sum(I0)
#         VARXAS1 = sumj[ VARIj/[sum(I)]^2 + VARI0j/[sum(I0)]^2 ]
#     XRF:
#         f' = 1
#         dXAS1/dIj  = 1 / sum(I0)
#         dXAS1/dI0j = -sum(I)/[sum(I0)]^2
#         VARXAS1 = sumj[ VARIj/[sum(I)]^2 + VARI0j/[sum(I0)]^2 ] . [sum(I)]^2/[sum(I0)]^2
#
# XAS2:
#     XAS2 = f(sum(I/I0)) = f(nrep.f^-1(m))
#          = nrep.m         (XRF)
#          = m - ln(nrep)   (Transmission)
#     VARXAS = sumj[ (dXAS/dIj)^2.VARIj + (dXAS/dI0j)^2.VARI0j ]
#
#     dXAS2/dIj  = f'/I0j
#     dXAS2/dI0j = -f'.Ij/T0j^2
#
#     Transmission:
#         f' = -1/sum(I/I0)
#         dXAS2/dIj  = -1/I0j/sum(I/I0)
#         dXAS2/dI0j = Ij/I0j^2/sum(I/I0)
#         VARXAS2 = sumj[ VARIj/I0j^2 + VARI0j.Ij^2/I0j^4] / [sum(I/I0)]^2
#
#     XRF:
#         f' = 1
#         dXAS2/dIj  = 1/I0j
#         dXAS2/dI0j = -Ij/T0j^2
#         VARXAS2 = sumj[ VARIj/I0j^2 + VARI0j.Ij^2/I0j^4 ]
#
# XAS3:
#     XAS3 = sum(f(I/I0)) = nrep.m
#     VARXAS = sumj[ (dXAS/dIj)^2.VARIj + (dXAS/dI0j)^2.VARI0j ]
#
#     dXAS3/dIj  = f'/I0j
#     dXAS3/dI0j = -f'.Ij/I0j^2
#
#     Transmission:
#         f' = -I0j/Ij
#         dXAS3/dIj  = -1/Ij
#         dXAS3/dI0j = 1/Ij0
#         VARXAS3 = sumj[ VARIj/Ij^2 + VARI0j/Ij0^2 ]
#
#     XRF:
#         f' = 1
#         dXAS3/dIj  = 1/I0j
#         dXAS3/dI0j = -Ij/I0j^2
#         VARXAS3 = sumj[ (VARIj/Ij^2 + VARI0j/I0j^2)  . Ij^2/Ij0^2 ]
#

#execfile("initcctbx.py")

import os, sys
sys.path.insert(1,'/data/id21/inhouse/wout/dev/SpectroCrunch')

import numpy as np
import matplotlib.pyplot as plt

from spectrocrunch.materials.compoundfromcif import compoundfromcif as compoundf
from spectrocrunch.materials.mixture import mixture as mixturef
from spectrocrunch.materials.types import fractionType

from uncertainties import ufloat
from uncertainties.umath import exp as uexp
from uncertainties.umath import exp as ulog

import xraylib

import warnings
warnings.filterwarnings("ignore")

def genxas(xrf=False,fine=False,refresh=True):
    # Material
    compound1 = compoundf("cinnabar",name='cinnabar')
    compound2 = compoundf("gypsum",name='gypsum')
    mixture = mixturef([compound1,compound2],[0.5,0.5],fractionType.weight)
    thickness = 10 # micron

    # Energies
    n = 120
    energy = np.linspace(2.46,2.52,n)

    # XRF
    if xrf:
        shells = [xraylib.K_SHELL]
        fluolines = [xraylib.__dict__[s] for s in xraylib.__dict__.keys() if s.endswith("_LINE") and s.startswith("K")]
        mixture.markabsorber("S",shells=shells,fluolines=fluolines)
        mu = mixture.partial_mass_abs_coeff(energy,decomposed=False,fine=fine,refresh=refresh)
        collection = 0.01
        m = mu*mixture.density()*thickness*1e-4*collection
    else:
        mu = mixture.mass_att_coeff(energy,decomposed=False,fine=fine,refresh=refresh)
        m = mu*mixture.density()*thickness*1e-4

    if xrf:
        fim = m
    else:
        fim = np.exp(-m)

    nrepeats = 10

    flux = np.full((nrepeats,n),1e5)
    data = np.round(flux * np.tile(fim,(nrepeats,1)))

    return energy,data,flux,m

def xasnorm(data,flux,xrf=False,normtype=1):
    nrepeats,n = data.shape

    # Assume Poisson
    fluxwnoise = np.random.poisson(lam=flux, size=None).astype(float)
    datawnoise = np.random.poisson(lam=data, size=None).astype(float)
    varI = data*1.5+100.
    varI0 = flux*1.5

    uind = -1
    udata = [ufloat(v1,np.sqrt(v2)) for v1,v2 in zip(data[:,uind],varI[:,uind])]
    uflux = [ufloat(v1,np.sqrt(v2)) for v1,v2 in zip(flux[:,uind],varI0[:,uind])]

    if normtype==1:
        # XAS1 = f(sum(I)/sum(I0))
        # VARXAS1 = sumj[ VARIj/[sum(I)]^2 + VARI0j/[sum(I0)]^2 ]
        # VARXAS1 = sumj[ VARIj/[sum(I)]^2 + VARI0j/[sum(I0)]^2 ] . [sum(I)]^2/[sum(I0)]^2
        print("XAS1 = sum(f(I/I0)")
        xas = data.sum(axis=0)/flux.sum(axis=0)
        xasnoise = datawnoise.sum(axis=0)/fluxwnoise.sum(axis=0)
        uxas = sum(udata)/sum(uflux)

        sumIsq = data.sum(axis=0)**2
        sumI0sq = flux.sum(axis=0)**2
        var = varI.sum(axis=0)/sumIsq + varI0.sum(axis=0)/sumI0sq

        if xrf:
            var *= sumIsq/sumI0sq
        else:
            xas = -np.log(xas)
            xasnoise = -np.log(xasnoise)
            uxas = -ulog(uxas)

    elif normtype==2:
        # XAS2 = f(sum(I/I0))
        # VARXAS2 = sumj[ VARIj/I0j^2 + VARI0j.Ij^2/I0j^4] / [sum(I/I0)]^2
        # VARXAS2 = sumj[ VARIj/I0j^2 + VARI0j.Ij^2/T0j^4 ]
        print("XAS2 = sum(f(I/I0)")
        xas = (data/flux).sum(axis=0)
        xasnoise = (datawnoise/fluxwnoise).sum(axis=0)
        uxas = sum([v1/v2 for v1,v2 in zip(udata,uflux)])

        Isq = data**2
        I0sq = flux**2
        var = ((varI + varI0*Isq/I0sq)/I0sq).sum(axis=0)

        if xrf:
            xas /= nrepeats
            xasnoise /= nrepeats
            uxas /= nrepeats
            var /= nrepeats**2
        else:
            xas = -np.log(xas) + np.log(nrepeats)
            xasnoise = -np.log(xasnoise) + np.log(nrepeats)
            uxas = -ulog(uxas) + np.log(nrepeats)
            var /= (data/flux).sum(axis=0)**2

    elif normtype==3:
        # XAS3 = sum(f(I/I0)) 
        # VARXAS3 = sumj[ VARIj/Ij^2 + VARI0j/Ij0^2 ]
        # VARXAS3 = sumj[ (VARIj/Ij^2 + VARI0j/I0j^2)  . Ij^2/Ij0^2 ]
        print("XAS3 = sum(f(I/I0)")
        xas = data/flux
        xasnoise = datawnoise/fluxwnoise
        uxas = [v1/v2 for v1,v2 in zip(udata,uflux)]

        Isq = data**2
        I0sq = flux**2
        var = (varI/Isq + varI0/I0sq)

        if xrf:
            var *= Isq/I0sq
        else:
            xas = -np.log(xas)
            xasnoise = -np.log(xasnoise)
            uxas = [-ulog(v) for v in uxas]

        xas = xas.sum(axis=0)
        xas /= nrepeats
        xasnoise = xasnoise.sum(axis=0)
        xasnoise /= nrepeats
        uxas = sum(uxas)
        uxas /= nrepeats
        var = var.sum(axis=0)
        var /= nrepeats**2
    else:
        raise ValueError("Unknown normalization type.")

    print("{}+/-{}".format(xas[uind],np.sqrt(var[uind])))
    print(uxas)
    return xas,var,xasnoise

if __name__ == '__main__':

    xrf = True

    energy,data,flux,att = genxas(xrf=xrf)

    for normtype in range(1,4):
        xas,var,xasnoise = xasnorm(data,flux,xrf=xrf,normtype=normtype)

        fig1 = plt.figure(1)
        ax = fig1.add_subplot(111)
        ax.plot(energy,xasnoise,label="Type {}".format(normtype))
        ax.set_xlabel("Energy (keV)")
        ax.set_ylabel("Mu.rho.d")

        fig2 = plt.figure(2)
        ax = fig2.add_subplot(111)
        ax.plot(energy,np.sqrt(var)/xas*100,label="Type {}".format(normtype))
        ax.set_xlabel("Energy (keV)")
        ax.set_ylabel("Noise/signal (%)")

    fig1 = plt.figure(1)
    ax = fig1.add_subplot(111)
    ax.plot(energy,att,label="Theory")
    plt.legend()

    fig2 = plt.figure(2)
    plt.legend()

    #plt.show()


