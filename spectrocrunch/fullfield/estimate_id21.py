# -*- coding: utf-8 -*-
#
#   Copyright (C) 2018 European Synchrotron Radiation Facility, Grenoble, France
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

import os
from ..io import ascii
from ..utils import timing
from ..utils import units
from ..materials.element import Element
from ..data import axis
from ..simulation import calcnoise
from ..materials.types import fraction

import scipy.optimize
import numpy as np
import matplotlib.pyplot as plt

class FFsetup(calcnoise.id21_ffsetup):

    def __init__(self,sample,layerindex=0,compoundname=None):
        super(FFsetup,self).__init__(sample)
        self.layerindex = layerindex
        self.compoundname = compoundname

    def set_wcompound(self,wcompound):
        x = wcompound/100.
        if x==1:
            raise ValueError('Stay below 100%')
            
        fracs = self.sample[self.layerindex].massfractions()
        s = sum([v for k,v in fracs.items() if k!=self.compoundname])
        m = (1-x)/s
        fracs = {k:v*m for k,v in fracs.items()}
        fracs[self.compoundname] = x
        
        self.sample[self.layerindex].change_fractions(fracs,fraction.mass)

    def get_wcompound(self):
        return self.sample[self.layerindex].massfractions()[self.compoundname]*100
    
    def set_vcompound(self,vcompound):
        x = vcompound/100.
        if x==1:
            raise ValueError('Stay below 100%')
        
        fracs = self.sample[self.layerindex].volumefractions()
        s = sum([v for k,v in fracs.items() if k!=self.compoundname])
        m = (1-x)/s
        fracs = {k:v*m for k,v in fracs.items()}
        fracs[self.compoundname] = x

        self.sample[self.layerindex].change_fractions(fracs,fraction.volume)

    def get_vcompound(self):
        return self.sample[self.layerindex].volumefractions()[self.compoundname]*100
        
    def set_tcompound(self,tcompound):
        self.set_vcompound(tcompound/self.get_layerthickness()*100)
        
    def get_tcompound(self):
        return self.get_vcompound()*self.get_layerthickness()/100

    def set_layerthickness(self,layerthickness):   
        self.sample[self.layerindex].thickness = layerthickness*1e-4

    def get_layerthickness(self):
        return self.sample[self.layerindex].thickness*1e4
        
    def optimize_thickness(self,I0,**kwargs):
        def costfunc(layerthickness):
            self.set_layerthickness(layerthickness[0])
            c = self.costfunc(I0,energy,**kwargs)
            return c
        
        guess = self.get_layerthickness()
        result = scipy.optimize.least_squares(costfunc, guess, gtol=1e-015, ftol=1e-015)
        return result.x[0],result.success

    def optimize_wcompound(self,I0,energy,**kwargs):
        def costfunc(wcompound):
            self.set_wcompound(wcompound[0])
            c = self.costfunc(I0,energy,**kwargs)
            return c
        
        guess = self.get_wcompound()
        result = scipy.optimize.least_squares(costfunc, guess, bounds=([0,100]), gtol=1e-015, ftol=1e-015)
        return result.x[0],result.success
        
    def optimize_thickness_plot(self,I0,energy,**kwargs):
        thickness = self.get_layerthickness()
        
        t = np.linspace(max(thickness-100,0),thickness+100,50)
        r = np.zeros(len(t))
        for i,layerthickness in enumerate(t):
            self.set_layerthickness(layerthickness)
            r[i] = self.costfunc(I0,energy,**kwargs)

        self.set_layerthickness(thickness)

        plt.plot(t,1/r,label="{} wt%".format(self.get_wcompound()))
        plt.xlabel("thickness ($\mu$m)")
        plt.ylabel("Jump-to-noise")
        
    def optimize_wcompound_plot(self,I0,energy,**kwargs):
        w = self.get_wcompound()
    
        t = np.linspace(0,99,100)
        r = np.zeros(len(t))
        for i,wcompound in enumerate(t):
            self.set_wcompound(wcompound)
            r[i] = self.costfunc(I0,energy,**kwargs)

        self.set_wcompound(w)
        
        plt.plot(t,1/r,label="{} $\mu$m".format(self.get_layerthickness()))
        plt.xlabel("{} (wt%)".format(self.compoundname))
        plt.ylabel("Jump-to-noise")
    
    def optimize_tcompound_plot(self,I0,energy,**kwargs):
        thickness = self.get_tcompound()
    
        tlayer = self.get_layerthickness()
    
        t = np.linspace(0,tlayer*0.99,100)
        r = np.zeros(len(t))
        for i,tcompound in enumerate(t):
            self.set_tcompound(tcompound)
            r[i] = self.costfunc(I0,energy,**kwargs)

        self.set_tcompound(thickness)
        
        plt.plot(t,1/r,label="{} $\mu$m".format(self.get_layerthickness()))
        plt.xlabel("{} thickness ($\mu$m)".format(self.compoundname))
        plt.ylabel("Jump-to-noise")
    
    def optimize(self,I0,energy,**kwargs):
        def costfunc(p):
            self.set_wcompound(p[0])
            self.set_layerthickness(p[1])
            return self.costfunc(I0,energy,**kwargs)
        
        guess = (self.get_wcompound(),self.get_layerthickness())
        
        result = scipy.optimize.least_squares(costfunc, guess, bounds=([0,0],[100,1e6]), gtol=1e-015)
        return result.x,result.success
    
    @property
    def edge_energies(self):
        return units.Quantity(self.element.edge_energies(),'keV').to('keV')

    def energybeforeafter(self,mi=-1,ma=+1):
        limits = units.Quantity([mi,ma],'eV') + self.edge_energies[0]
        return limits.to('keV').magnitude
    
    def JNR(self,flux,**kwargs):
        return 1/self.costfunc(flux,self.energybeforeafter(**kwargs),**self.simulkwargs)[0]
    
    @property
    def energies(self):
        return self.axis.magnitude
    
    def setconfig(self,element=None,edge=None,limits=None,stepsizes=None,ndark=None,
                  timeperenergy=None,flatfraction=None,tframe_data=None,tframe_flat=None):
        self.element = Element(element)
        self.element.markabsorber(shells = edge)
        limits = units.Quantity(limits,'eV') + self.edge_energies[0]
        limits.ito('keV')
        stepsizes = units.Quantity(stepsizes,'eV').to('keV')
        self.axis = axis.AxisSegments(limits,[10]*(len(limits)-1))
        self.axis.stepsizes = stepsizes
        self.timeinfo = [{'timeperenergy':tm,'flatfraction':frac,'tframe_data':tfd,'tframe_flat':tff}
                          for tm,frac,tfd,tff in zip(timeperenergy,flatfraction,tframe_data,tframe_flat)]
        self.ndark = ndark

    @property
    def simulkwargs(self):
        tframe_data = np.mean([d['tframe_data'] for d in self.timeinfo])
        tframe_flat = np.mean([d['tframe_flat'] for d in self.timeinfo])
        ndata = np.mean([d['timeperenergy']/d['tframe_data']*(1-d['flatfraction']) for d in self.timeinfo])
        nflat = np.mean([d['timeperenergy']/d['tframe_flat']*d['flatfraction'] for d in self.timeinfo])
        ndata = int(round(ndata))
        nflat = int(round(nflat))
        return {"tframe_data":tframe_data,"nframe_data":ndata,\
                    "tframe_flat":tframe_flat,"nframe_flat":nflat,\
                    "nframe_dark":self.ndark}
    
    @property
    def simul_timeperenergy(self):
        kwargs = self.simulkwargs
        return kwargs['nframe_data']*kwargs['tframe_data'] + kwargs['nframe_flat']*kwargs['tframe_flat']
        
    @property
    def simul_flatfraction(self):
        kwargs = self.simulkwargs
        return kwargs['nframe_flat']*kwargs['tframe_flat']/self.simul_timeperenergy
        
    def genconfig(self,path,name):
        limits = self.axis.limits
        nsteps = self.axis.nsteps
        stepsizes = self.axis.stepsizes

        fmth = '{:<10}: {:d}'
        fmtf = '  #{:d} {:<6}: {:.04f}'
        fmtd = '  #{:d} {:<6}: {:d}'
        fmtr = 'Range: {:~.04f} - {:~.04f} in {:d} x {:~.03f}'
        fmtt = ' {}: {}'
        with ascii.Writer(os.path.join(path,name+'.cfg')) as fcfg:
            with ascii.Writer(os.path.join(path,name+'.txt')) as ftxt:
                
                fcfg.write(fmth.format('dark',1))
                fcfg.write(fmth.format('dark nb im',self.ndark))
                fcfg.write(fmth.format('refafter',1))
                fcfg.write(fmth.format('nb. zone',len(nsteps)))
                
                tmdatatot = 0
                tmflattot = 0
                tmoverheadtot = 0
                for i,(n,stepsize,info) in enumerate(zip(nsteps,stepsizes,self.timeinfo)):
                    
                    ndata,nflat = self.getnframes(info['timeperenergy'],info['tframe_data'],
                                                  info['tframe_flat'],info['flatfraction'])
                    tmdata = ndata*info['tframe_data']
                    tmflat = nflat*info['tframe_flat']
                    tmoverhead = 6.50305 + 0.0131498*(ndata+nflat)
                    
                    tmdatatot += tmdata*n
                    tmflattot += tmflat*n
                    tmoverheadtot += tmoverhead*n
                    
                    fcfg.write(fmtf.format(i,'start',limits[i].to('keV').magnitude))
                    fcfg.write(fmtf.format(i,'end',limits[i+1].to('keV').magnitude))
                    fcfg.write(fmtd.format(i,'nbp',n))
                    fcfg.write(fmtf.format(i,'time',info['tframe_data']))
                    fcfg.write(fmtd.format(i,'frame',ndata))
                    fcfg.write(fmtf.format(i,'time',info['tframe_flat']))
                    fcfg.write(fmtd.format(i,'frame',nflat))

                    ftxt.write(fmtr.format(limits[i],limits[i+1],n,stepsize.to('eV')))
                    ftxt.write(fmtt.format('Total time',timing.strseconds((tmdata+tmflat+tmoverhead)*n)))
                    ftxt.write(fmtt.format(' data',timing.strseconds(tmdata*n)))
                    ftxt.write(fmtt.format(' flat',timing.strseconds(tmflat*n)))
                    ftxt.write(fmtt.format(' overhead',timing.strseconds(tmoverhead*n)))
                
                ftxt.write('\nNumber of energies: {}'.format(len(self.axis)))
                ftxt.write('Total time: {}'.format(timing.strseconds(tmdatatot+tmflattot+tmoverheadtot)))
                ftxt.write(fmtt.format('data',timing.strseconds(tmdatatot)))
                ftxt.write(fmtt.format('flat',timing.strseconds(tmflattot)))
                ftxt.write(fmtt.format('overhead',timing.strseconds(tmoverheadtot)))
