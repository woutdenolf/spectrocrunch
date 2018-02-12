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

from . import base
from ..materials import compoundfromname
from ..materials import element
from ..common import constants
from ..common.classfactory import with_metaclass
from ..common import instance

import scipy.special
import numpy as np

class Detector(with_metaclass(base.PointSourceCentric)):

    def __init__(self,mcazero=None,mcagain=None,mcanoise=None,mcafano=None,\
                    shape_fixedarearatios=None,shape_pymca=None,shape_conversionenergy=None,**kwargs):
        
        super(Detector,self).__init__(**kwargs)
        
        self.mcazero = mcazero # keV
        self.mcagain = mcagain # keV/bin
        self.mcanoise = mcanoise # FWHM in keV
        self.mcafano = mcafano
        
        self.shape_conversionenergy = shape_conversionenergy
        self.fixedarearatios = shape_fixedarearatios is not None
        if self.fixedarearatios:
            self.tailbroadening= shape_fixedarearatios["tailbroadening"]
            self.fractions = (shape_fixedarearatios["tailfraction"],shape_fixedarearatios["stepfraction"])
        else:
            if shape_pymca is not None:
                self.ratios = (shape_pymca["tailarea_ratio"],shape_pymca["stepheight_ratio"])
                self.tailslope_ratio = shape_pymca["tailslope_ratio"]

    @property
    def tailbroadening(self):
        return self._tailbroadening
    
    @tailbroadening.setter
    def tailbroadening(self,value):
        self._tailbroadening = value
        if self.fixedarearatios and self.shape_conversionenergy is not None:
            self.tailslope_ratio = self._conv_tailslope_ratio()
    
    @property
    def fractions(self):
        wtail = self._tailfraction
        wstep = self._stepfraction
        wdet = 1-wtail-wstep
        return wdet,wtail,wstep
    
    @fractions.setter
    def fractions(self,value):
        wtail,wstep = value
        if wtail+wstep>1:
            raise RuntimeError("Tail and step fractions must be <= 1")
        self._tailfraction = wtail
        self._stepfraction = wstep
        
        if self.fixedarearatios and self.shape_conversionenergy is not None:
            self.ratios = self._conv_ratios()
    
    @property
    def tailslope_ratio(self):
        return self._tailslope_ratio
    
    @tailslope_ratio.setter
    def tailslope_ratio(self,value):
        self._tailslope_ratio = value
        
        if not self.fixedarearatios and self.shape_conversionenergy is not None:
            self.tailbroadening = self._conv_tailbroadening()
    
    @property
    def ratios(self):
        return self._tailarea_ratio,self._stepheight_ratio

    @ratios.setter
    def ratios(self,value):
        self._tailarea_ratio,self._stepheight_ratio = value

        if not self.fixedarearatios and self.shape_conversionenergy is not None:
            self.fractions = self._conv_fractions()

    def _conv_tailslope_ratio(self):
        gvar = self.gaussianVAR(self.shape_conversionenergy)
        return self.tailbroadening*np.sqrt(gvar)
    
    def _conv_tailbroadening(self):
        gvar = self.gaussianVAR(self.shape_conversionenergy)
        return self.tailslope_ratio/np.sqrt(gvar)

    def _conv_ratios(self):
        # rtail = wtail/wdet
        # rstep = wstep/wdet*gnorm/energy
        
        wdet,wtail,wstep = self.fractions
        if wdet==0:
            wdet = 1
            
        gvar = self.gaussianVAR(self.shape_conversionenergy)
        gnorm = np.sqrt(2*np.pi*gvar)
        c = gnorm/self.shape_conversionenergy
        rtail = wtail/wdet
        rstep = wstep/wdet*c
        return rtail,rstep

    def _conv_fractions(self):
        # rtail = wtail/(1-wtail-wstep)
        # rstep = wstep/(1-wtail-wstep)*c
        #
        # rtail = rtail    * wstep + (rtail+1)* wtail
        # rstep = (rstep+c)* wstep + rstep    * wtail
        

        rtail,rstep = self.ratios
        
        gvar = self.gaussianVAR(self.shape_conversionenergy)
        gnorm = np.sqrt(2*np.pi*gvar)
        return self._calc_fractions(rtail,rstep,self.shape_conversionenergy,gnorm)

    def _calc_fractions(self,rtail,rstep,energy,gnorm,nopeak=False):
        c = gnorm/energy
        if rtail==0 and rstep==0:
            wtail = 0
            wstep = 0
        elif rtail==0:
            wtail = 0
            wstep = rstep/(rstep+c)
        elif rstep==0:
            wtail = rtail/(rtail+1.)
            wstep = 0
        else:
            ctail = (rtail+1.)/rtail
            if rtail > rstep:
                wtail = c/( (rstep+c)*ctail - rstep)
                wstep = 1 - ctail*wtail
            else:
                cstep = (rstep+c)/rstep
                wstep = (1 - ctail) / (1 - cstep/ctail)
                wtail = 1 - cstep*wstep

        return wtail,wstep
        
    def __str__(self):
        if self.fixedarearatios:
            shape = " Tail broadening = {}\n"\
                    " Peak fraction = {}\n"\
                    " Tail fraction = {}\n"\
                    " Step fraction = {}"\
                    .format(self.tailbroadening,*self.fractions)
        else:
            shape = " Tail slope ratio = {}\n"\
                    " Tail area ratio = {}\n"\
                    " Step height ratio = {}"\
                    .format(self.tailslope_ratio,*self.ratios)

        return "XRF detector:\n{}\n"\
                "MCA:\n zero = {} eV\n"\
                " gain = {} eV\n"\
                " noise = {} eV (FWHM)\n"\
                " fano = {}\n{}"\
                .format(super(Detector,self).__str__(),\
                self.mcazero*1000,self.mcagain*1000,\
                self.mcanoise*1000,self.mcafano,shape)

    def FWHMtoVAR(self,FWHM):
        return FWHM**2/(8*np.log(2))

    def VARtoFWHM(self,var):
        return np.sqrt(var*(8*np.log(2)))
        
    def gaussianVAR(self,energy):
        """Gaussian variance (keV^2)
        """
        return self.FWHMtoVAR(self.mcanoise) + self.mcafano*self.ehole*1e-3*energy

    def gaussianFWHM(self,energy):
        """Gaussian FWHM (keV)
        """
        return self.VARtoFWHM(self.gaussianVAR(energy))

    def voigtFWHM(self,energy,linewidth=0):
        fG = self.gaussianFWHM(energy) # detection FWHM
        fL = linewidth # transition FWHM
        return 0.5346*fL + np.sqrt(0.2166*fL**2+fG**2)
    
    def voigtVAR(self,energy):
        return self.FWHMtoVAR(self.voigtFWHM(energy))

    def lineprofile(self,x,u,linewidth=0):
        """
        Args:
            x(num|array): energies (keV, nx)
            u(num|array): peak energies (keV, nu)
            linewidth(Optional(num|array)): natural line widths
            
        Returns:
            array: nx x nu
        """
        # [1] R. van Grieken, A. Markowicz: Handbook of X-ray Spectrometry, CRC Press, Boca Raton (2001)
        #     peak, tail and step have fixed area ratio's
        # [2] PyMca
        #     peak and tail have a fixed area ratio
        #     peak and step have a fixed height ratio
        
        x = instance.asarray(x)[:,np.newaxis]
        u = instance.asarray(u)[np.newaxis,:] 
        diff = x-u

        gvar = self.gaussianVAR(u)
        a = 2*gvar # 2.sigma^2
        gnorm = np.sqrt(np.pi*a) # sqrt(2.pi.sigma^2)

        if self.fixedarearatios:
            tailbroadening = self.tailbroadening
            wdet,wtail,wstep = self.fractions
            bdet = wdet>0
            btail = tailbroadening>0 and wtail>0
            bstep = wstep>0
        else:
            tailslope_ratio = self.tailslope_ratio
            tailarea_ratio,stepheight_ratio = self.ratios
            wtail,wstep = self._calc_fractions(tailarea_ratio,stepheight_ratio,u,gnorm)
            wdet = 1-wtail-wstep
            bdet = wdet>0
            btail = tailslope_ratio>0 and tailarea_ratio>0
            bstep = stepheight_ratio>0
        
            #wdet = 1. # I think this is what PyMca does
        
        if bdet or bstep or btail:
            b = np.sqrt(a) # sqrt(2).sigma
            if bstep or btail:
                argstep = diff/b # (x-u)/(sqrt(2).sigma)

        # Detector response:
        #   Gaussian: H*exp(-(x-u)^2/(2.gvar))
        #   Lorentz:  2/(pi.W)/(1+4/W^2.(x-u)^2)    (W:FWHM)
        #             W/(2.pi)/(W^2/4+(x-u)^2)
        #             y/pi/(y^2+(x-u)^2)            (y:W/2)
        #   Voigt:    Re(w[(x-u+i.y)/(2.sqrt(gvar))])/sqrt(2.pi.gvar)
        #
        #   VanGrieken: H = wdet/sqrt(2.pi.gvar)
        #   Pymca: H = garea/sqrt(2.pi.gvar)
        #   => garea = wdet
        if bdet:
            height = wdet/gnorm
            W = instance.asarray(linewidth)[:,np.newaxis]
            if W.any():
                y = height*np.real(scipy.special.wofz((diff + 0.5j*W)/b))
            else:
                y = height*np.exp(-diff**2/a)
        else:
            y = np.zeros_like(diff)
             
        # Incomplete charge collection (step):
        #   Gaussian: H.erfc[(x-u)/sqrt(2.gvar)]
        #
        #   VanGrieken: H = wstep/(2.u)
        #   Pymca: H = stepheight_ratio*garea/(2.gnorm)
        #   => stepheight_ratio = wstep/wdet*gnorm/u
        if bstep:
            # wstep/u = wdet/gnorm*stepheight_ratio/gnorm
            if self.fixedarearatios:
                semiheight = wstep/(2.*u)
            else:
                semiheight = stepheight_ratio*wdet/(2.*gnorm)
            y += semiheight * scipy.special.erfc(argstep)

        # Incomplete charge collection (tail):
        #   Gaussian: exp[(x-u)/(tb.sqrt(gvar))+1/(2.tb**2)].erfc[(x-u)/sqrt(2.gvar)+1/(tb.sqrt(2))]/(2.tb.sqrt(gvar))
        #             tr = tb.sqrt(gvar)
        #             H.exp[(x-u)/tr+gvar/(2.tr**2)].erfc[(x-u)/sqrt(2.gvar)+sqrt(gvar)/(sqrt(2).tr)]
        #
        #   VanGrieken: H = wtail/(2.tr)
        #   Pymca: H = garea*tailarea_ratio/(2.tr)
        #   => tailarea_ratio = wtail/wdet
        if btail:
            if self.fixedarearatios:
                tailslope_ratio = tailbroadening*np.sqrt(gvar)
                m = wtail/(2.*tailslope_ratio)
            else:
                m = tailarea_ratio*wdet/(2.*tailslope_ratio)
            
            with np.errstate(over="ignore"):
                mexp = m*np.exp(diff/tailslope_ratio+gvar/(2.*tailslope_ratio**2))
                ind = np.isinf(mexp)
                if ind.any():
                    mexp[ind] = 0 # lim_x->inf exp(x)*erfc(x) = 0
            y += mexp*scipy.special.erfc(argstep+b/(2.*tailslope_ratio))
            
        return y
    
    def addtopymca(self,setup,cfg): 
        super(Detector,self).addtopymca(setup,cfg)

        mcazero = self.mcazero
        mcagain = self.mcagain
        mcanoise = self.mcanoise
        # Compensate fano for difference in Ehole with pymca's 3.85 eV (ClassMcaTheory)
        mcafano = self.mcafano * self.ehole / 3.85
        cfg["detector"]["zero"] = mcazero
        cfg["detector"]["gain"] = mcagain
        cfg["detector"]["noise"] = mcanoise
        cfg["detector"]["fano"] = mcafano

        # No one-to-one correspondance
        energy = np.max(setup.energy)

        rtail,rstep = self.ratios
        bhypermet = True
        bstail = rtail>0
        bltail = False
        bstep = rstep>0
        
        if bstail:
            cfg["peakshape"]["st_arearatio"] = rtail
            cfg["peakshape"]["st_sloperatio"] = self.tailslope_ratio
        elif bltail:
            cfg["peakshape"]["lt_arearatio"] = rtail
            cfg["peakshape"]["lt_sloperatio"] = self.tailslope_ratio
            
        if bstep:
            cfg["peakshape"]["step_heightratio"] = rstep
        
        cfg["fit"]["hypermetflag"] = 1*bhypermet | 2*bstail | 4*bltail | 8*bstep

        xmin = int(round((setup.emin-mcazero)/mcagain))
        xmax = int(round((setup.emax-mcazero)/mcagain))
        cfg["fit"]["xmin"] = xmin
        cfg["fit"]["xmax"] = xmax
        cfg["fit"]["use_limit"] = 1
 
    def loadfrompymca(self,setup,cfg):
        super(Detector,self).loadfrompymca(setup,cfg)
        
        self.activearea = cfg["concentrations"]["area"]
        self.distance = cfg["concentrations"]["distance"]
        self.mcazero = cfg["detector"]["zero"]
        self.mcagain = cfg["detector"]["gain"]
        self.mcanoise = cfg["detector"]["noise"]
        # Compensate fano for difference in Ehole with pymca's 3.85 eV (ClassMcaTheory)
        self.mcafano = cfg["detector"]["fano"] * 3.85 / self.ehole
        
        # No one-to-one correspondance
        energy = np.max(setup.energy)

        b = cfg["fit"]["hypermetflag"]
        bstail = (b&2)==2
        bltail = (b&4)==4
        bstep = (b&8)==8

        if bstail:
            rtail = cfg["peakshape"]["st_arearatio"]
            self.tailslope_ratio = cfg["peakshape"]["st_sloperatio"]
        elif bltail:
            rtail = cfg["peakshape"]["lt_arearatio"]
            self.tailslope_ratio = cfg["peakshape"]["lt_sloperatio"]
        else:
            rtail = 0
            self.tailslope_ratio = 0
            
        if bstep:
            rstep = cfg["peakshape"]["step_heightratio"]
        else:
            rstep = 0
            
        self.ratios = (rtail,rstep)

        setup.emin = cfg["fit"]["xmin"]*self.mcagain+self.mcazero
        setup.emax = cfg["fit"]["xmax"]*self.mcagain+self.mcazero

class sn3102(Detector):
    aliases = ["XFlash5100"]
    
    def __init__(self,**kwargs):
        ultralene = compoundfromname.compoundfromname("ultralene")
        moxtek = compoundfromname.compoundfromname("moxtek ap3.3")
        
        attenuators = {}
        attenuators["FoilDetector"] = {"material":ultralene,"thickness":4e-4} # cm
        attenuators["WindowDetector"] = {"material":moxtek,"thickness":380e-4} # AP3.3 Moxtek
        attenuators["Detector"] = {"material":element.Element('Si'),"thickness":450e-4}
        kwargs["attenuators"] = attenuators
        
        kwargs["activearea"] = 80e-2*0.75 # cm^2
        
        kwargs["mcazero"] = 0. # keV
        kwargs["mcagain"] = 5e-3 # keV
        kwargs["mcanoise"] = 0.1 # keV
        kwargs["mcafano"] = 0.114
        
        #kwargs["shape_fixedarearatios"] = {"tailbroadening":0.5,"tailfraction":0.05,"stepfraction":0.005}
        kwargs["shape_pymca"] = {"stepheight_ratio":0.001,"tailarea_ratio":0.05,"tailslope_ratio":0.5}
        
        kwargs["ehole"] = constants.eholepair_si()

        super(sn3102,self).__init__(**kwargs)
        
              
class Leia(Detector):
    aliases = ["SGX80"]
    
    def __init__(self,**kwargs):
        ultralene = compoundfromname.compoundfromname("ultralene")
        
        attenuators = {}
        attenuators["FoilDetector"] = {"material":ultralene,"thickness":4e-4} #cm
        attenuators["WindowDetector"] = {"material":element.Element('Be'),"thickness":25e-4}
        attenuators["Detector"] = {"material":element.Element('Si'),"thickness":450e-4}
        kwargs["attenuators"] = attenuators
        
        kwargs["activearea"] = 80e-2 # cm^2
        
        kwargs["mcazero"] = 0. # keV
        kwargs["mcagain"] = 5e-3 # keV
        kwargs["mcanoise"] = 50e-3 # keV
        kwargs["mcafano"] = 0.19
        
        #kwargs["shape_fixedarearatios"] = {"tailbroadening":0.5,"tailfraction":0.05,"stepfraction":0.005}
        kwargs["shape_pymca"] = {"stepheight_ratio":0.001,"tailarea_ratio":0.05,"tailslope_ratio":0.5}

        kwargs["ehole"] = constants.eholepair_si()

        super(Leia,self).__init__(**kwargs)


class BB8(Detector):
    aliases = ["SGX50"]
    
    def __init__(self,**kwargs):
        ultralene = compoundfromname.compoundfromname("ultralene")
        
        attenuators = {}
        attenuators["FoilDetector"] = {"material":ultralene,"thickness":4e-4} #cm
        attenuators["WindowDetector"] = {"material":element.Element('Be'),"thickness":12.5e-4}
        attenuators["Detector"] = {"material":element.Element('Si'),"thickness":450e-4}
        kwargs["attenuators"] = attenuators
        
        kwargs["activearea"] = 50e-2 # cm^2
        
        kwargs["mcazero"] = 0. # keV
        kwargs["mcagain"] = 5e-3 # keV
        kwargs["mcanoise"] = 0.1 # keV
        kwargs["mcafano"] = 0.114
        
        #kwargs["shape_fixedarearatios"] = {"tailbroadening":0.5,"tailfraction":0.05,"stepfraction":0.005}
        kwargs["shape_pymca"] = {"stepheight_ratio":0.001,"tailarea_ratio":0.05,"tailslope_ratio":0.5}
        
        kwargs["ehole"] = constants.eholepair_si()
        
        super(BB8,self).__init__(**kwargs)


class DR40(Detector):
    aliases = ["VITUSH80"]
    
    def __init__(self,**kwargs):
        ultralene = compoundfromname.compoundfromname("ultralene")
        
        attenuators = {}
        attenuators["FoilDetector"] = {"material":ultralene,"thickness":4e-4} #cm
        attenuators["WindowDetector"] = {"material":element.Element('Be'),"thickness":25e-4}
        attenuators["Detector"] = {"material":element.Element('Si'),"thickness":450e-4}
        kwargs["attenuators"] = attenuators
        
        kwargs["activearea"] = 80e-2 # cm^2
        
        kwargs["mcazero"] = 0. # keV
        kwargs["mcagain"] = 5e-3 # keV
        kwargs["mcanoise"] = 0.1 # keV
        kwargs["mcafano"] = 0.114
        
        #kwargs["shape_fixedarearatios"] = {"tailbroadening":0.5,"tailfraction":0.05,"stepfraction":0.005}
        kwargs["shape_pymca"] = {"stepheight_ratio":0.001,"tailarea_ratio":0.05,"tailslope_ratio":0.5}
        
        kwargs["ehole"] = constants.eholepair_si()
        
        super(DR40,self).__init__(**kwargs)
        
factory = Detector.factory

