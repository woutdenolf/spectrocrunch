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

from ..simulation.classfactory import with_metaclass
from ..simulation import noisepropagation

class Detector(with_metaclass()):

    def __init__(self, etoDU=None, qe=None, DUoffset=None, darkcurrent=None, readoutnoise=None, geometry=None,**kwargs):
        """
        Args:
            etoDU(num): number of DU per electron (DU/e)
            qe(num): detector quantum efficiency (e/ph)
            DUoffset(num): pixel intensity offset (DU)
            darkcurrent(num): dark current (e/sec)
            readoutnoise(num): readout noise (e)
            geometry()
        """
        self.etoDU = etoDU
        self.qe = qe
        self.DUoffset = DUoffset
        self.darkcurrent = darkcurrent
        self.readoutnoise = readoutnoise
        self.geometry = geometry
        super(Detector,self).__init__(**kwargs)
        
    def __str__(self):
        return "qe = {} e-/DU\n conv = {} e-/DU\n off = {} DU\n dark = {} e-/sec\n noise = {} e-\n{}".format(self.qeself.etoDU,self.DUoffset,self.darkcurrent,self.readoutnoise,str(self.geometry))
    
    def propagate(self,N,visspectrum,tframe=None,nframe=None,forward=True,poissonapprox=False):
        """Error propagation of a number of photons.
               
        Args:
            N(unumpy.uarray): incomming number of photons within a tframe timespan with uncertainties
            visspectrum(spectrum): visible light spectrum
            tframe(num|numpy.array): time per frame (sec)
            nframe(num|numpy.array): number of frames (sec)

        Returns:
            numpy.array:
        """

        qe = self.qe(visspectrum) 

        N,qe = self.propagate_broadcast(N,qe)

        if noisepropagation.israndomvariable(N):
            darkcurrent = noisepropagation.poisson(self.darkcurrent)
            readoutnoise = noisepropagation.randomvariable(0,self.readoutnoise)
            
            process = noisepropagation.poisson(qe)
            if forward:
                Nout = noisepropagation.compound(N,process,forward=forward)
                
                if poissonapprox:
                    ENout = noisepropagation.E(Nout)
                    Nout = noisepropagation.randomvariable(ENout,ENout**0.5)
                    
                Nout = (Nout + darkcurrent*tframe + readoutnoise)*self.etoDU + self.DUoffset
                Nout = noisepropagation.repeat(nframe,Nout,forward=forward)
            else:
                Nout = noisepropagation.repeat(nframe,N,forward=forward)
                Nout = noisepropagation.reverse_add(Nout,self.DUoffset)
                Nout = noisepropagation.reverse_mult(Nout,self.etoDU)
                Nout = noisepropagation.reverse_add(Nout,darkcurrent*tframe)
                Nout = noisepropagation.reverse_add(Nout,readoutnoise)
                
                if poissonapprox:
                    ENout = noisepropagation.E(Nout)
                    Nout = noisepropagation.randomvariable(ENout,ENout**0.5)
                
                Nout = noisepropagation.compound(Nout,process,forward=forward)
        else:
            if forward:
                Nout = nframe* ((N*qe + self.darkcurrent*tframe)*self.etoDU + self.DUoffset)
            else:
                Nout = ((N/nframe - self.DUoffset)/self.etoDU - self.darkcurrent*tframe)/qe
        
        return Nout
        
class pcoedge55(Detector):
    """
    PCO Edge 5.5
    """
    aliases = ["PCO Edge 5.5"]

    def __init__(self,geometry=None):
        qe = lambda visspectrum: 0.7
        super(pcoedge55, self).__init__(etoDU=65536/30000., qe=qe, DUoffset=95.5, darkcurrent=7.4, readoutnoise=0.95, geometry=geometry)

class frelon2k16(Detector):
    """
    Frelon 2K 16
    """
    aliases = ["Frelon 2K 16"]

    def __init__(self,geometry=None):
        # https://doi.org/10.1063/1.2783112
        # https://doi.org/10.1107/S0909049506000550
        qe = lambda visspectrum: 0.7
        super(frelon2k16, self).__init__(etoDU=65536/320000., qe=qe, DUoffset=100., darkcurrent=1., readoutnoise=19., geometry=geometry)

class frelon2k14(Detector):
    """
    Frelon 2K 14
    """
    aliases = ["Frelon 2K 14"]

    def __init__(self,geometry=None):
        # https://doi.org/10.1063/1.2783112
        # https://doi.org/10.1107/S0909049506000550
        qe = lambda energy: [0.7]*energy.size
        super(frelon2k14, self).__init__(etoDU=16384/320000., qe=qe, DUoffset=100., darkcurrent=1., readoutnoise=24., geometry=geometry)

classes = Detector.clsregistry
aliases = Detector.aliasregistry
factory = Detector.factory


