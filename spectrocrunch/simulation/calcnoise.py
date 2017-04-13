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

from . import detectors
from . import scintillators
from . import lenses

def id21_ffnoise(I0,energy,sample,tframe,nframe):
    """ID21 fullfield noise propagation

    Args:
        I0 (num or array-like): incomming flux (ph/sec)
        energy (num or array-like): associated energy (keV)
        sample (spectrocrunch.simulation.materials.Material): sample
        tframe(num): time per frame (sec)
        nframe(num): number of frames (sec)

    Returns:
        uncertainties.unumpy.uarray: detector signal in ADU
    """

    N = I0*tframe

    oscint = scintillators.Scintillator.factory("LSO ID21",10)
    olens = lenses.Lens.factory("mitutoyoid21_10x")
    odet = detectors.AreaDetector.factory("pcoedge55")

    N = sample.propagate(N,energy)
    N = oscint.propagate(N,energy)
    N = olens.propagate(N,energy,nrefrac=oscint.get_nrefrac())
    N = odet.propagate(N,energy,tframe=tframe,nframe=nframe)

    return N

