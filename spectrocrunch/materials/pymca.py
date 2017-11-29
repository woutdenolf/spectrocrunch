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

import compoundfromformula
import compoundfromlist
import compoundfromname
import mixture
import element
import types
from ..common.classfactory import with_metaclass
from ..common import units

import itertools


############################################################

class Geom(with_metaclass(object)):

    def __init__(self,anglein=None,angleout=None,detectorposition=None):
        self.anglein = float(anglein) # deg
        self.angleout = float(angleout) # deg
        self.detectorposition = float(detectorposition) # mm

class sdd60(Geom):

    def __init__(self,**kwargs):
        super(sdd60,self).__init__(anglein=62,angleout=49,**kwargs)

    @property
    def distance(self):
        return self.detectorposition + 60.5 # blc10516 (April 2017)

class sdd90(Geom):

    def __init__(self,**kwargs):
        super(sdd90,self).__init__(anglein=62,angleout=28,**kwargs)

    @property
    def distance(self):
        return self.detectorposition + 85.5 # blc10516 (April 2017)

############################################################

class Detector(with_metaclass(object)):

    def __init__(self,mcazero=None,mcagain=None,activearea=None,geometry=None,**kwargs):
        self.mcazero = mcazero # eV
        self.mcagain = mcagain # eV/bin
        self.activearea = float(activearea) # mm^2
        self.geom = Geom.factory(geometry,**kwargs)

    @property
    def distance(self):
        return self.geom.distance
        
    @property
    def solidangle(self):
        distance = self.distance
        r2 = self.activearea/np.pi # squared disk radius
        return 2.*np.pi*(1.-(distance/np.sqrt(r2+distance**2.)))

    def distancefromsolidangle(self,soldiangle):
        solidanglefrac = soldiangle/(4*np.pi)
        if solidanglefrac>=0.5:
            raise ValueError("Solid angle must be < 2.pi")
        r2 = self.activearea/np.pi # squared disk radius
        self.geom.distance = np.sqrt(r2)*(0.5-solidanglefrac)/np.sqrt((1-solidanglefrac)*solidanglefrac)

    def addtoconfig(self,config,energy): 
        config["concentrations"]["area"] = self.activearea/100. # mm^2 -> cm^2
        config["concentrations"]["distance"] = self.distance/10. # mm -> cm
        
        mcazero = self.mcazero/1000. # eV -> keV
        mcagain = self.mcagain/1000.
        config["detector"]["zero"] = mcazero
        config["detector"]["gain"] = mcagain
        
        emin = 0.9
        emax = energy+0.5
        xmin = int((emin-mcazero)/mcagain)
        xmax = int((emax-mcazero)/mcagain)
        config["fit"]["xmin"] = xmin
        config["fit"]["xmax"] = xmax
        config["fit"]["use_limit"] = 1

class Leia(Detector):

    def __init__(self,**kwargs):
        super(Leia,self).__init__(activearea=80,**kwargs)
        
    def addtoconfig(self,config,energy):
        super(Leia,self).addtoconfig(config,energy)
        config["attenuators"]["kapton"] = [1, 'Mylar', 1.4, 4e-4, 1.0]
        config["attenuators"]["window"] = [1, 'Be1', element.Element('Be').density, 25e-4, 1.0]
        config["attenuators"]["Detector"] = [1, 'Si1', element.Element('Si').density, 0.5, 1.0]

class BB8(Detector):

    def __init__(self,**kwargs):
        super(BB8,self).__init__(activearea=80,**kwargs)
        
    def addtoconfig(self,config,energy):
        super(BB8,self).addtoconfig(config,energy)
        config["attenuators"]["kapton"] = [1, 'Mylar', 1.4, 4e-4, 1.0]
        config["attenuators"]["window"] = [1, 'Be1', element.Element('Be').density, 12.5e-4, 1.0]
        config["attenuators"]["Detector"] = [1, 'Si1', element.Element('Si').density, 0.5, 1.0]

        
############################################################

class Sample(with_metaclass(object)):
    
    def __init__(self,layers=None,thickness=None,names=None,environment=None,detector=None,**kwargs):
        self.layers = layers
        self.thickness = thickness
        self.names = names
        self.environment = environment
        self.detector = Detector.factory(detector,**kwargs)

    def __iter__(self):
        return itertools.izip(self.layers,self.thickness,self.names)

    def addmaterial(self,config,mat,thickness,name=None):
        print ""
        print name
        for k,v in mat.weightfractions().items():
            print k,v
        k,v = mat.pymcaformat(thickness=thickness,name=name)
        config["materials"][k] = v
        return k
        
    def defmatrix(self,config,name):
        anglein = self.detector.geom.anglein
        angleout = self.detector.geom.angleout 
        if name=="MULTILAYER":
            density,thickness = 0.,0.
        else:
            v = config["materials"][name]
            density,thickness = v["Density"], v["Thickness"]
        config["attenuators"]["Matrix"] = [1, name, density, thickness, anglein, angleout, 0, anglein+angleout]
    
    def deflayer(self,config,index,name):
        v = config["materials"][name]
        layer = "Layer{}".format(index)
        config["multilayer"][layer] = [1,name,v["Density"],v["Thickness"]]

    def addshells(self,config,elements,energy):
        emax = energy
        emin = 0.9
        
        for e in elements:
            shells = element.Shell.pymcafactory((e,emin,emax))
            if shells:
                config["peaks"][str(e)] = shells
                
    def addtoconfig(self,config,energy):
        self.detector.addtoconfig(config,energy)

        if len(self.layers)==1:
            layer,t,name = iter(self).next()
            name = self.addmaterial(config,layer,t,name=name)
            self.defmatrix(config,name)
        else:
            for index,(layer,t,name) in enumerate(self):
                name = self.addmaterial(config,layer,t,name=name)
                self.deflayer(config,index,name)
 
            self.defmatrix(config,'MULTILAYER')
        
        for l in self.layers:
            self.addshells(config,l.elements,energy)

        if self.environment is not None:
            self.addshells(config,self.environment,energy)

def axo(name,elements,ad,windowthickness,filmthickness):
    # When the filmthickness is not known exactly, there is no other
    # option than assume all elements are in the substrate. In this case
    # we would want to switch off absorption corrections.
    if filmthickness is None:
        w = compoundfromname.compoundfromname("silicon nitride")
        arealdensity = w.arealdensity()

        arealdensity.update(dict(zip(elements,ad)))
        
        totalarealdensity = sum(arealdensity.values())
        massfractions = {e:ad/totalarealdensity for e,ad in arealdensity.items()}
        
        layer1 = compoundfromlist.CompoundFromList(massfractions.keys(),massfractions.values(),types.fractionType.weight,w.density)
        
        return [layer1],[windowthickness*1e-7],[name]
    
    else:
        elements = [compoundfromlist.CompoundFromList([e],[1],types.fractionType.mole,0,name=e) for e in elements]
        layer1 = mixture.Mixture(elements,ad,types.fractionType.weight)

        layer2 = compoundfromname.compoundfromname("silicon nitride")
        
        return [layer1,layer2],[filmthickness*1e-7,windowthickness*1e-7],[name,"silicon nitride"]

class AXOID21_1(Sample):
    aliases = ["RF7-200-S2371-03"]
    
    def __init__(self,**kwargs):
        name = "RF7-200-S2371-03"
        elements = ["Pb","La","Pd","Mo","Cu","Fe","Ca"]
        ad = [7.7,9,1.9,0.9,2.4,4,11.4]
        windowthickness = 200 # nm
        filmthickness = None # nm
        layers,thickness,names = axo(name,elements,ad,windowthickness,filmthickness)
        
        super(AXOID21_1,self).__init__(layers=layers,thickness=thickness,names=names,**kwargs)


class AXOID21_2(Sample):
    aliases = ["RF8-200-S2454"]
    
    def __init__(self,**kwargs):
        name = "RF8-200-S2454"
        elements = ["Pb","La","Pd","Mo","Cu","Fe","Ca"]
        ad = [6.3,7.6,2.3,0.7,2.6,4.1,25.1]
        windowthickness = 200 # nm
        filmthickness = None # nm
        layers,thickness,names = axo(name,elements,ad,windowthickness,filmthickness)
        
        super(AXOID21_1,self).__init__(layers=layers,thickness=thickness,names=names,**kwargs)

class AXOID16b_1(Sample):
    aliases = ["RF8-200-S2453"]
    
    def __init__(self,**kwargs):
        name = "RF8-200-S2453"
        elements = ["Pb","La","Pd","Mo","Cu","Fe","Ca"]
        ad = [5.9,10.3,1.2,.7,2.4,3.9,20.3]
        windowthickness = 200 # nm
        filmthickness = None # nm
        layers,thickness,names = axo(name,elements,ad,windowthickness,filmthickness)
        
        super(AXOID16b_1,self).__init__(layers=layers,thickness=thickness,names=names,**kwargs)
        
############################################################

factory = Sample.factory


class PymcaConfig(object):

    def __init__(self,sample=None,energy=None,flux=None,time=None,fluolinefilter=False):
        self.sample = sample
        self.energy = energy # keV
        self.flux = flux # ph/s
        self.time = time # s
        self.config = {}
    
    @property
    def _energy(self):
        return units.magnitude(self.energy,"keV")
    
    @property
    def _time(self):
        return units.magnitude(self.time,"s")
    
    @property
    def _flux(self):
        return units.magnitude(self.flux,"hertz")
        
    def beam(self,config):
        energy = self._energy
        if False:
            config["fit"]["energy"] = [energy, 100*energy]
            config["fit"]["energyflag"] = [1,1-self.fluolinefilter]
            config["fit"]["energyscatter"] = [1,0]
            config["fit"]["energyweight"] = [1e+100,1e-05]
        else:
            config["fit"]["energy"] = [energy]
            config["fit"]["energyflag"] = [1]
            config["fit"]["energyscatter"] = [1]
            config["fit"]["energyweight"] = [1.0]
        config["fit"]["scatterflag"] = 1
        
        config["concentrations"]["flux"] = self._flux
        config["concentrations"]["time"] = self._time
        
    def background(self,config):
        config["fit"]["stripflag"] = 1
        config["fit"]["stripalgorithm"] = 1
        config["fit"]["snipwidth"] = 100
    
    def peakshape(self,config):
        config["fit"]["hypermetflag"] = 3 # shorttail
        config["fit"]["hypermetflag"] = 1

    def fill(self,config):
        self.sample.addtoconfig(config,self._energy)
        self.beam(config)
        self.background(config)
        self.peakshape(config)


