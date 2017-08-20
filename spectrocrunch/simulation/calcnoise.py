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

from __future__ import division

from . import areadetectors
from . import scintillators
from . import lenses
from .noisepropagation import randomvariable

import numpy as np

def id21_ffnoise(I0,energy,sample,\
                    tframe_data=1,nframe_data=1,\
                    tframe_flat=1,nframe_flat=1,\
                    nframe_dark=1,\
                    scint="LSO",lens="x10"):
    """ID21 fullfield noise propagation

    Args:
        I0 (num or array-like): incomming flux (ph/sec)
        energy (num or array-like): associated energy (keV)
        sample (spectrocrunch.simulation.materials.Material): sample
        tframe_data(num): time per frame (sec)
        nframe_data(num): number of data frames
        tframe_flat(num): time per frame (sec)
        nframe_flat(num): number of flat frames
        nframe_dark(num): number of dark frames
        scint(Optional(str)):
        lens(Optional(str)):

    Returns:
        4-tuple(uncertainties.unumpy.uarray): detector signal in ADU (with and w.o. sample, with and w.o. exposure)
    """

    if scint=="LSO":
        oscint = scintillators.factory("LSO ID21",thickness=10)
    else:
        oscint = scintillators.factory("GGG ID21",thickness=13)

    if lens=="x10":
        olens = lenses.factory("mitutoyoid21_10x")
    else:
        olens = lenses.factory("mitutoyoid21_10x")

    odet = areadetectors.factory("pcoedge55")

    # With sample
    Ni = I0*tframe_data
    Ni = randomvariable(Ni,np.sqrt(Ni))
    Ni = sample.propagate(Ni,energy)
    Ni = oscint.propagate(Ni,energy)
    Ni = olens.propagate(Ni,energy,nrefrac=oscint.get_nrefrac())
    Ni = odet.propagate(Ni,energy,tframe=tframe_data,nframe=nframe_data)
    
    # Without sample
    N0 = I0*tframe_flat
    N0 = randomvariable(N0,np.sqrt(N0))
    N0 = oscint.propagate(N0,energy)
    N0 = olens.propagate(N0,energy,nrefrac=oscint.get_nrefrac())
    N0 = odet.propagate(N0,energy,tframe=tframe_flat,nframe=nframe_flat)

    # Without beam
    Di = odet.propagate(randomvariable(0,0),energy,tframe=tframe_data,nframe=nframe_dark)
    Di *= nframe_data/nframe_dark
    D0 = odet.propagate(randomvariable(0,0),energy,tframe=tframe_flat,nframe=nframe_dark)
    D0 *= nframe_flat/nframe_dark

    return Ni,N0,Di,D0

def id21_transmissionnoise(I0,energy,time,iodetgain,idetgain,iodet="1"):
    """ID21 transmission noise propagation
    """
    pass

def id21_xrfnoise(I0,energy,time):
    """ID21 XRF noise propagation
    """
    pass


