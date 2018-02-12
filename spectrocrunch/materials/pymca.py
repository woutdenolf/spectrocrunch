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

from ..common import units
from ..common import instance
from . import mixture
from . import element
from .. import xraylib

from . import compoundfromformula

import fisx
import contextlib
import numpy as np
import re
import scipy.interpolate
import copy
from PyMca5.PyMcaPhysics.xrf import ClassMcaTheory
from PyMca5.PyMcaPhysics.xrf import ConcentrationsTool
from PyMca5.PyMcaIO import ConfigDict

class PymcaHandle(object):

    def __init__(self,sample=None,emin=None,emax=None,\
                energy=None,weights=None,scatter=None,\
                flux=None,time=None,escape=1,ninteractions=1):
        self.sample = sample
        
        self.energy = units.magnitude(energy,"keV")
        self._emin = emin
        self._emax = emax
        self.weights = weights
        self.scatter = scatter
        self.escape = escape
        self.ninteractions = ninteractions
        
        self.flux = units.magnitude(flux,"hertz")
        self.time = units.magnitude(time,"s")
        
        self.mcafit = ClassMcaTheory.McaTheory()
        self.ctool = ConcentrationsTool.ConcentrationsTool()
    
    def __str__(self):
        s = zip(instance.asarray(self.energy),instance.asarray(self.weights),instance.asarray(self.scatter))
        s = '\n '.join("{} keV: {} % (Scatter: {})".format(k,v*100,sc) for k,v,sc in s)
        return "Flux = {} ph/s\nTime = {} s\nSource lines:\n {}\n{}\n{}".format(self.flux,self.time,s,self.sample,self.sample.geometry)
    
    @property
    def emax(self):
        if self._emax is None:
            # include elastic scattering peaks
            p = np.max(self.energy)
            s = np.sqrt(self.sample.geometry.detector.gaussianVAR(p))
            return p+3*s
        else:
            return self._emax

    @emax.setter
    def emax(self,value):
        self._emax = value
        
    @property
    def emin(self):
        if self._emin is None:
            # include Al-K lines
            p = xraylib.LineEnergy(13,xraylib.KL3_LINE)
            s = np.sqrt(self.sample.geometry.detector.gaussianVAR(p))
            return p-3*s
        else:
            return self._emin

    @emin.setter
    def emin(self,value):
        self._emin = value

    def addtopymca_lstsq(self,cfg):
        cfg["fit"]["maxiter"] = 500
      
    def addtopymca_beam(self,cfg):
        # TODO: move to source?
        cfg["fit"]["energy"] = instance.asarray(self.energy).tolist()
        cfg["fit"]["energyweight"] = instance.asarray(self.weights).tolist()
        cfg["fit"]["energyscatter"] = instance.asarray(self.scatter,dtype=int).tolist()
        cfg["fit"]["energyflag"] = np.ones_like(cfg["fit"]["energy"],dtype=int).tolist()
        cfg["fit"]["scatterflag"] = 1
        
        cfg["concentrations"]["flux"] = self.flux
        cfg["concentrations"]["time"] = self.time
    
    def loadfrompymca_beam(self,cfg):
        ind = np.asarray(cfg["fit"]["energyflag"],dtype=bool).tolist()
    
        self.energy = instance.asarray(cfg["fit"]["energy"])
        self.energy = self.energy[ind]
        self.weights = instance.asarray(cfg["fit"]["energyweight"])
        self.weights = self.weights[ind]
        self.scatter = instance.asarray(cfg["fit"]["energyscatter"])
        self.scatter = self.scatter[ind]
        self.flux = cfg["concentrations"]["flux"]
        self.time = cfg["concentrations"]["time"]

    def addtopymca_background(self,cfg):
        rtail,rstep = self.sample.geometry.detector.ratios
        cfg["fit"]["stripflag"] = int(rstep<=0)
        cfg["fit"]["stripalgorithm"] = 1
        cfg["fit"]["snipwidth"] = 100
    
    def addtopymca_other(self,cfg):
        cfg["fit"]["escapeflag"] = 1
        cfg["concentrations"]["usemultilayersecondary"] = self.ninteractions-1
    
    def loadfrompymca_other(self,cfg):
        self.ninteractions = cfg["concentrations"]["usemultilayersecondary"]+1
        
    def addtopymca_material(self,cfg,material,defaultthickness=1e-4):
        matname,v = material.topymca(defaultthickness=defaultthickness)
        if material.nelements>1:
            cfg["materials"][matname] = v
        return matname
    
    def loadfrompymca_material(self,cfg,matname,density):
        if matname in cfg["materials"]:
            material = mixture.Mixture.frompymca(cfg["materials"][matname])
        else:
            pattern = "^(?P<element>[A-Z][a-z]?)1$"
            m = re.match(pattern,matname)
            if m:
                m = m.groupdict()
                material = element.Element(m["element"])
            else:
                material = compoundfromformula.CompoundFromFormula(matname,density)
        return material
        
    def setdata(self,y):
        x = np.arange(len(y), dtype=np.float32)
        self.mcafit.setData(x,y)

    def loadfrompymca(self,filename=None,config=None):
        if filename is not None:
            fconfig = ConfigDict.ConfigDict()
            fconfig.read(filename)
            self.mcafit.configure(fconfig)
            
        if config is not None:
            self.mcafit.configure(config)
        
        config = self.mcafit.getConfiguration()
        self.loadfrompymca_beam(config)
        self.loadfrompymca_other(config)
        self.sample.loadfrompymca(self,config)
    
    def addtopymca(self,fresh=False,onlygeometry=False):
        if fresh:
            config = self.mcafit.getStartingConfiguration()
            config["attenuators"] = {}
            try:
                self.mcafit.useFisxEscape(True)
            except:
                pass
            self.mcafit.disableOptimizedLinearFit()
        else:
            config = self.mcafit.getConfiguration()
        
        if onlygeometry:
            self.sample.geometry.addtopymca(self,config)
        else:
            self.addtopymca_beam(config)
            self.addtopymca_background(config)
            self.addtopymca_lstsq(config)
            self.addtopymca_other(config)
            self.sample.addtopymca(self,config)
        
        self.mcafit.configure(config)

    def processfitresult(self,digestedresult,originalconcentrations=False):
        ctoolcfg = self.ctool.configure()
        ctoolcfg.update(digestedresult['config']['concentrations'])
        return self.ctool.processFitResult(config=ctoolcfg,
                                                        fitresult={"result":digestedresult},
                                                        elementsfrommatrix=originalconcentrations,
                                                        fluorates = self.mcafit._fluoRates,
                                                        addinfo=True)
    
    def concentrationsfromfitresult(self,digestedresult,out=None):
        # Mass fractions:
        #
        # For a group (e.g Fe-K) and i loops over the layers:
        #  grouparea = flux.time.solidangle/(4.pi).sum_i[massfrac_i.grouprate_i]
        #
        # When Fe is present in only one layer j:
        #  massfrac_j = grouparea/(flux.time.solidangle/(4.pi).grouprate_i)
        #
        # When Fe present in multiple layers:
        #  grouparea = flux.time.solidangle/(4.pi).massfrac_avg.sum_i[b_i.grouprate_i]
        # where b_i is 1 or 0 depending on whether it is present in the matrix definition
        # When the element is not in the matrix definition b_i=1 for all layers.
        #
        #  massfrac_avg = grouparea/(flux.time.solidangle/(4.pi).sum_i[b_i.grouprate_i])
        
        if out is None:
            out = {}
        
        conresult,addinfo = self.processfitresult(digestedresult,originalconcentrations=False)
        
        out["massfractions"] = {k.replace(' ','-'):v for k,v in conresult["mass fraction"].items()}
        
        if "layerlist" in conresult:
            out["lmassfractions"] = [{k.replace(' ','-'):v for k,v in conresult[k]["mass fraction"].items()} for k in conresult["layerlist"]]
        if len(out["lmassfractions"])==0:
            out["lmassfractions"] = [out["massfractions"]]

        # Group rates:
        #   grouprate = solidangle/(4.pi).sum_i[massfrac_i.grouprate_i]
        
        if True:
            out["rates"] = {}
            safrac = addinfo["SolidAngle"]
            #assert(self.sample.geometry.detector.solidangle/(4*np.pi)==addinfo["SolidAngle"])
            #assert(self.flux*self.time==addinfo["I0"])
            for group in out["massfractions"]:
                element,shell = group.split("-")
                grouprate = 0.
                for i in range(1,self.sample.nlayers+1):
                    if element in self.mcafit._fluoRates[i]:
                        massfrac_l = self.mcafit._fluoRates[i][element]["mass fraction"]
                        grouprate_l = self.mcafit._fluoRates[i][element]["rates"]["{} xrays".format(shell)]
                        grouprate += massfrac_l*grouprate_l
                grouprate *= safrac
                out["rates"][group] = grouprate

        out["fitareas"] = {k.replace(' ','-'):v for k,v in conresult["fitarea"].items()}
        
        out["fitrates"] = {k:v/addinfo["I0"] for k,v in out["fitareas"].items()}
        
        return out
    
    def spectrafromfitresult(self,digestedresult,out=None):
        conresult,addinfo = self.processfitresult(digestedresult,originalconcentrations=True)

        # Group rates:
        if False:
            # TODO: set rates of elements not in matrix to zero
            out["rates"] = {}
            for group in out["area"]:
                out["rates"][group.replace(' ','-')] = conresult['area'][group]/addinfo["I0"]
        
        # Matrix spectrum
        nparams = len(digestedresult['parameters'])
        ngroups = len(digestedresult['groups'])
        nglobal = nparams - ngroups
        parameters = list(digestedresult['fittedpar'])
        for i,group in enumerate(digestedresult['groups']):
            if group in conresult['area']:
                # Replace fitted with theoretical peak area
                parameters[nglobal+i] = conresult['area'][group]
            else:
                # Scattering peaks don't have a theoretical peak area
                parameters[nglobal+i] = 0.
        ymatrix = self.mcafit.mcatheory(parameters,digestedresult['xdata'])

        yback = digestedresult["continuum"]
        
        gen = lambda spectrum:scipy.interpolate.interp1d(digestedresult["energy"],spectrum,\
                    kind="nearest",bounds_error=False,fill_value=(spectrum[0],spectrum[-1]))

        if out is None:
            out = {}
        out["energy"] = digestedresult["energy"]
        out["channels"] = digestedresult["xdata"]
        out["y"] = digestedresult["ydata"]
        out["yfit"] = digestedresult["yfit"]
        out["yback"] = yback
        out["gen"] = gen
        out["ypileup"] = digestedresult["pileup"]
        out["ymatrix"] = ymatrix+yback
                
        return out
    
    def configfromfitresult(self,digestedresult):
        nparams = len(digestedresult['parameters'])
        ngroups = len(digestedresult['groups'])
        nglobal = nparams - ngroups
        
        parammap = {"Zero":("detector","zero"),\
                    "Gain":("detector","gain"),\
                    "Noise":("detector","noise"),\
                    "Fano":("detector","fano"),\
                    "Sum":("detector","sum"),\
                    "ST AreaR":("peakshape","st_arearatio"),\
                    "ST SlopeR":("peakshape","st_sloperatio"),\
                    "LT AreaR":("peakshape","lt_arearatio"),\
                    "LT SlopeR":("peakshape","lt_sloperatio"),\
                    "STEP HeightR":("peakshape","step_heightratio"),\
                    "Eta Factor":("peakshape","eta_factor")}
        
        config = copy.deepcopy(digestedresult['config'])
        
        for i in range(nglobal):
            param = digestedresult['parameters'][i]
            if param in parammap:
                key = parammap[param]
                config[key[0]][key[1]] = digestedresult['fittedpar'][i]
        
        return config
        
    def fit(self):
        # Fit
        self.mcafit.estimate()
        fitresult,digestedresult = self.mcafit.startfit(digest=1)

        # Load parameters from fit
        self.loadfrompymca(config=self.configfromfitresult(digestedresult))
        
        # Parse result
        result = {}
        self.concentrationsfromfitresult(digestedresult,out=result)
        self.spectrafromfitresult(digestedresult,out=result)
        
        return result

    def xrayspectrum(self,**kwargs):
        return self.sample.xrayspectrum(self.energy,emin=self.emin,emax=self.emax,\
                    weights=self.weights,ninteractions=self.ninteractions,**kwargs)
        
    def xraygrouprates(self,**kwargs):
        spectrum = self.xrayspectrum(**kwargs)
        groups = spectrum.linegroups(convert=True)
        pattern = "^(?P<element>[A-Z][a-z]?)-(?P<group>[A-Z]).?$"
        groups2 = {}
        for k in groups:
            m = re.match(pattern,k)
            if m:
                m = m.groupdict()
                g = "{}-{}".format(m["element"],m["group"])
                A = sum(groups[k].values())
                if g in groups2:
                    groups2[g] += A
                else:
                    groups2[g] = A
    
        return groups2


class FisxConfig():
    
    FISXMATERIALS = None
    
    def init(self):
        if self.FISXMATERIALS is None:
            self.FISXMATERIALS = fisx.Elements()
            self.FISXMATERIALS.initializeAsPyMca()
        
    @contextlib.contextmanager
    def init_ctx(self):
        self.init()
        yield
        
    def setDetector(self,setup,detector):
        with self.init_ctx():
            detector.addtofisx(setup,self.FISXMATERIALS)

    def addtofisx_material(self,material):
        if material.nelements>1:
            with self.init_ctx():
                self.FISXMATERIALS.addMaterial(material.tofisx(),errorOnReplace=False)
        return material.pymcaname

