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

from .simul import with_simulmetaclass

from . import noisepropagation

import numpy as np

class AreaDetector(with_simulmetaclass()):
    """
    Class representing an area detector
    """

    def __init__(self, etoDU=None, qe=None, DUoffset=None, darkcurrent=None, readoutnoise=None):
        """
        Args:
            etoDU(num): number of DU per electron (DU/e)
            qe(num): detector quantum efficiency (e/ph)
            DUoffset(num): pixel intensity offset (DU)
            darkcurrent(num): dark current (e/sec)
            readoutnoise(num): readout noise (e)
        """
        self.required(etoDU,"etoDU")
        self.required(qe,"qe")
        self.required(DUoffset,"DUoffset")
        self.required(darkcurrent,"darkcurrent")
        self.required(readoutnoise,"readoutnoise")
        
        self.etoDU = float(etoDU)
        self.qe = qe
        
        self.DUoffset = float(DUoffset)
        
        self.darkcurrent = noisepropagation.poisson(darkcurrent)
        
        # Read-out noise
        # https://spie.org/samples/PM170.pdf
        self.readoutnoise = noisepropagation.randomvariable(0,readoutnoise)

    def propagate(self,N,visspectrum,tframe=None,nframe=None,withnoise=True,forward=True):
        """Error propagation of a number of photons.
               
        Args:
            N(unumpy.uarray): incomming number of photons within a tframe timespan with uncertainties
            visspectrum(spectrum): visible light spectrum
            tframe(num|numpy.array): time per frame (sec)
            nframe(num|numpy.array): number of frames (sec)

        Returns:
            numpy.array:
        """

        self.defined(self.propagate,tframe,"frame exposure time")
        self.defined(self.propagate,nframe,"number of frames")

        qe = self.qe(visspectrum) 
        
        N,qe = self.propagate_broadcast(N,qe)
        
        if withnoise:
            process = noisepropagation.poisson(qe)

            if forward:
                Nout = noisepropagation.compound(N,process,forward=forward)
                Nout = (Nout + self.darkcurrent*tframe)*self.etoDU + self.DUoffset
                Nout = noisepropagation.repeat(nframe,Nout,forward=forward)
            else:
                Nout = noisepropagation.repeat(nframe,N,forward=forward)
                Nout = noisepropagation.reverse_add(Nout,self.DUoffset)
                Nout = noisepropagation.reverse_mult(Nout,self.etoDU )
                Nout = noisepropagation.reverse_add(Nout,self.darkcurrent*tframe )
                Nout = noisepropagation.compound(Nout,process,forward=forward)
        else:
            if forward:
                Nout = nframe* ((N*qe + noisepropagation.E(self.darkcurrent)*tframe)*self.etoDU + self.DUoffset)
            else:
                Nout = ((N/nframe - self.DUoffset)/self.etoDU - noisepropagation.E(self.darkcurrent)*tframe)/qe
        
        return Nout

class pcoedge55(AreaDetector):
    """
    PCO Edge 5.5
    """
    aliases = ["PCO Edge 5.5"]

    def __init__(self):
        qe = lambda visspectrum: 0.7
        super(pcoedge55, self).__init__(etoDU=65536/30000., qe=qe, DUoffset=95.5, darkcurrent=7.4, readoutnoise=0.95)

class frelon2k16(AreaDetector):
    """
    Frelon 2K 16
    """
    aliases = ["Frelon 2K 16"]

    def __init__(self):
        # https://doi.org/10.1063/1.2783112
        # https://doi.org/10.1107/S0909049506000550
        qe = lambda visspectrum: 0.7
        super(pcoedge55, self).__init__(etoDU=65536/320000., qe=qe, DUoffset=100., darkcurrent=1., readoutnoise=19.)

class frelon2k14(AreaDetector):
    """
    Frelon 2K 14
    """
    aliases = ["Frelon 2K 14"]

    def __init__(self):
        # https://doi.org/10.1063/1.2783112
        # https://doi.org/10.1107/S0909049506000550
        qe = lambda energy: [0.7]*energy.size
        super(pcoedge55, self).__init__(etoDU=16384/320000., qe=qe, DUoffset=100., darkcurrent=1., readoutnoise=24.)
        
classes = AreaDetector.clsregistry
aliases = AreaDetector.aliasregistry
factory = AreaDetector.factory

