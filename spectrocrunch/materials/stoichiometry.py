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

import xraylib
import fdmnes

import numpy as np
import os
import tempfile

# compounds in mixture:
#
#  mass of compound i in a mixture
#  mi = Vi*rhoi(Vi) = vi*V*rhoi(Vi)
#
#  density of compound i in the mixture
#  rhoi(V) = mi/V = vi*rhoi(Vi)
#          = mi/V = wi*rho(V)
#
#  volume fraction as a function of weight fractions and densities
#  vi = Vi/total(Vj) = mi/rhoi(Vi) / total(mj/rhoj(Vj)) = wi/rhoi(Vi) / total(wj/rhoj(Vj))
#
#  weight fraction as a function of volume fractions and densities
#  wi = mi/total(mj) = vi*rhoi(Vi) / total(vj*rhoi(Vj))
#
#  mole fraction as a function of weight fractions and relative atomic weights
#  xi = ni/total(nj) = mi/MMi / total(mj/MMj) = wi/MMi / total(wj/MMj)
#
#  weight fraction as a function of mole fractions
#  wi = mi/total(mj) = ni*MMi / total(nj*MMj) = xi*MMi / total(xj*MMj)
#
#  density of the mixture as a function of volume fractions
#  rho(V) = total(mi)/V = total(vi*rhoi(Vi))
#
#  density of the mixture as a function of weight fractions
#  rho(V) = total(mi)/V = total(vi*rhoi(Vi)) = total(wi / total(wj/rhoj(Vj))) = 1/total(wj/rhoj(Vj))
#
#  mass attenuation coefficient of a mixture
#  mu = total(wi*mui)
#
#  linear attenuation coefficient of a mixture
#  muL = rho(V) * mu = total(wi*mui) / total(wj/rhoj(Vj))
#                    = total(rhoi(V)*mui)
#  
#  Molar mass of a compound
#  MM = M/n = total(mi)/total(mi/MMi) = total(wi)/total(wi/MMi) = 1/total(wi/MMi)

def frac_mole_to_weight(nfrac,MM):
    """
    Args:
        nfrac(np.array): mole fraction of each compound
        MM(np.array): molar mass of each compound
    """
    return(nfrac*MM/(nfrac*MM).sum())

def frac_weight_to_mole(wfrac,MM):
    """
    Args:
        wfrac(np.array): weight fraction of each compound
        MM(np.array): molar mass of each compound
    """
    return(wfrac/(MM * (wfrac/MM).sum()))

def molarmass_from_wfrac(wfrac,MM):
    """
    Args:
        wfrac(np.array): weight fraction of each compound
        MM(np.array): molar mass of each compound
    """
    return(1/((wfrac/MM).sum()))

def frac_weight_to_volume(wfrac,rho):
    """
    Args:
        wfrac(np.array): weight fraction of each compound
        rho(np.array): density of each compound
    """
    return(wfrac/(rho * (wfrac/rho).sum()))

def frac_volume_to_weight(vfrac,rho):
    """
    Args:
        vfrac(np.array): volume fraction of each compound
        rho(np.array): density of each compound
    """
    return(vfrac*rho / (vfrac*rho).sum())

def frac_volume_to_mole(vfrac,rho,MM):
    """
    Args:
        vfrac(np.array): volume fraction of each compound
        rho(np.array): density of each compound
        MM(np.array): molar mass of each compound
    """
    return(vfrac*rho/(MM * (vfrac*rho/MM).sum()))

def frac_mole_to_volume(nfrac,rho,MM):
    """
    Args:
        nfrac(np.array): mole fraction of each compound
        rho(np.array): density of each compound
        MM(np.array): molar mass of each compound
    """
    return(nfrac*MM/(rho * (nfrac*MM/rho).sum()))

def density_from_volumefrac(vfrac,rho):
    """
    Args:
        vfrac(np.array): volume fraction of each compound
        rho(np.array): density of each compound
    """
    return((vfrac*rho).sum())

def density_from_weightfrac(wfrac,rho):
    """
    Args:
        wfrac(np.array): weight fraction of each compound
        rho(np.array): density of each compound
    """
    return(1/(wfrac/rho).sum())

def density_from_molefrac(nfrac,rho,MM):
    """
    Args:
        nfrac(np.array): mole fraction of each compound
        rho(np.array): density of each compound
        MM(np.array): molar mass of each compound
    """
    return((nfrac*MM).sum()/(nfrac*MM/rho).sum())

