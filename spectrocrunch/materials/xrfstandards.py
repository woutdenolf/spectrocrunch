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

from . import compoundfromname
from . import compoundfromlist
from . import mixture
from . import element
from . import multilayer
from . import types
from . import pymca

from ..common.classfactory import with_metaclass

class Sample(with_metaclass(with_metaclass(pymca.PymcaAttenuators))):
    
    def __init__(self,material=None,environment=None,detector=None,**kwargs):
        if material is None:
            material = compoundfromname.compoundfromname("vacuum")
        self.material = material
        self.environment = environment
        self.detector = detector
        super(Sample,self).__init__(**kwargs)

    def __getattr__(self, attr):
        return getattr(self.detector,attr)
        
    def addtoconfig(self,config,energy):
        super(Sample,self).addtoconfig(config)
        self.detector.addtoconfig(config,energy)

        if isinstance(self.material,multilayer.Multilayer):
            for index,layer in enumerate(self.material):
                self.addlayer(config,index,layer)
                self.addshells(config,layer.elements,energy)
            self.addmatrix(config,'MULTILAYER')
        else:
            name = self.addmaterial(config,self.material)
            self.addshells(config,self.material.elements,energy)
            self.addmatrix(config,name)

        if self.environment is not None:
            self.addshells(config,self.environment,energy)

    def addmatrix(self,config,name):
        anglein = self.anglein
        angleout = self.angleout 
        if name=="MULTILAYER":
            density,thickness = 0.,0.
        else:
            v = config["materials"][name]
            density,thickness = v["Density"], v["Thickness"]
        config["attenuators"]["Matrix"] = [1, name, density, thickness, anglein, angleout, 0, anglein+angleout]
    
    def addlayer(self,config,index,layer):
        name = self.addmaterial(config,layer)
        l = "Layer{}".format(index)
        config["multilayer"][l] = [1,name,layer.density,layer.thickness]

    def addshells(self,config,elements,energy):
        emax = energy
        emin = 0.9
        
        for e in elements:
            shells = element.Shell.pymcafactory((e,emin,emax))
            if shells:
                config["peaks"][str(e)] = shells

    def loadfromconfig(self,config):
        super(Sample,self).loadfromconfig(config)
        self.detector.loadfromconfig(config,energy)
        
def axo(name,elements,ad,windowthickness,filmthickness):
    # When the filmthickness is not known exactly, there is no other
    # option than assume all elements are in the substrate. In this case
    # we would want to switch off absorption corrections.
    
    ultralene = compoundfromname.compoundfromname("ultralene")
        
    attenuators = {}
    attenuators["SampleCover"] = {"material":ultralene,"thickness":4e-4}
    attenuators["BeamFilter0"] = {"material":ultralene,"thickness":4e-4}

    if filmthickness is None:
        w = compoundfromname.compoundfromname("silicon nitride")
        arealdensity = w.arealdensity()

        arealdensity.update(dict(zip(elements,ad)))
        
        totalarealdensity = sum(arealdensity.values())
        massfractions = {e:ad/totalarealdensity for e,ad in arealdensity.items()}
        
        layer1 = compoundfromlist.CompoundFromList(massfractions.keys(),massfractions.values(),types.fractionType.weight,w.density,name=name)
        
        material = layer1
        #material = multilayer.Multilayer(layer1,windowthickness*1e-7)
    else:
        elements = [compoundfromlist.CompoundFromList([e],[1],types.fractionType.mole,0,name=e) for e in elements]
        layer1 = mixture.Mixture(elements,ad,types.fractionType.weight,name=name)

        layer2 = compoundfromname.compoundfromname("silicon nitride")
        
        material = multilayer.Multilayer([layer1,layer2],[filmthickness*1e-7,windowthickness*1e-7])
        
    return material,attenuators

class AXOID21_1(Sample):
    aliases = ["RF7-200-S2371-03"]
    
    def __init__(self,**kwargs):
        name = "RF7-200-S2371-03"
        elements = ["Pb","La","Pd","Mo","Cu","Fe","Ca"]
        ad = [7.7,9,1.9,0.9,2.4,4,11.4]
        windowthickness = 200 # nm
        filmthickness = None # nm
        material,attenuators = axo(name,elements,ad,windowthickness,filmthickness)

        super(AXOID21_1,self).__init__(material=material,attenuators=attenuators,**kwargs)


class AXOID21_2(Sample):
    aliases = ["RF8-200-S2454"]
    
    def __init__(self,**kwargs):
        name = "RF8-200-S2454"
        elements = ["Pb","La","Pd","Mo","Cu","Fe","Ca"]
        ad = [6.3,7.6,2.3,0.7,2.6,4.1,25.1]
        windowthickness = 200 # nm
        filmthickness = None # nm
        layers,thickness,attenuators = axo(name,elements,ad,windowthickness,filmthickness)
        
        super(AXOID21_1,self).__init__(material=material,attenuators=attenuators,**kwargs)

class AXOID16b_1(Sample):
    aliases = ["RF8-200-S2453"]
    
    def __init__(self,**kwargs):
        name = "RF8-200-S2453"
        elements = ["Pb","La","Pd","Mo","Cu","Fe","Ca"]
        ad = [5.9,10.3,1.2,.7,2.4,3.9,20.3]
        windowthickness = 200 # nm
        filmthickness = None # nm
        layers,thickness,attenuators = axo(name,elements,ad,windowthickness,filmthickness)
        
        super(AXOID16b_1,self).__init__(material=material,attenuators=attenuators,**kwargs)

factory = Sample.factory

