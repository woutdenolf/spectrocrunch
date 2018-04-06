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

class Standard(with_metaclass(multilayer.Multilayer)):
    
    def __init__(self,extra=None,**kwargs):
        if extra is None:
            self.extra = None
        else:
            self.extra = map(element.Element,extra)
        super(Standard,self).__init__(**kwargs)

    def addtopymca(self,setup,cfg):
        super(Standard,self).addtopymca(setup,cfg)
        if self.extra is not None:
            self.addtopymca_shells(setup,cfg,self.extra)


def axo(name,elements,arealdensity,windowthickness,filmthickness):
    """
    Args:
        name(str): name of standard
        elements(list(str)): element symbols
        arealdensity(list(num)): ng/mm^2
        windowthickness(num): substrate thickness (nm)
        filmthickness(num|None): deposited film thickness (nm)
        
    Returns:
        material(list(spectrocrunch.materials.compound|mixture)): layer composition
        thickness(list(num)): layer thickness in cm
        attenuators(list(list)): beam attenuators (before and after sample)
    """

    # When the filmthickness is not known exactly, there is no other
    # option than assume all elements are in the substrate. In this case
    # we would want to switch off absorption corrections.
    
    ultralene = compoundfromname.compoundfromname("ultralene")
    
    attenuators = [["SampleCover",ultralene,4e-4],\
                   ["BeamFilter0",ultralene,4e-4]]
    
    if filmthickness is None:
        arealdensity = dict(zip(elements,arealdensity))
        
        w = compoundfromname.compoundfromname("silicon nitride")
        arealdensity.update(w.arealdensity())
        
        totalarealdensity = sum(arealdensity.values())
        massfractions = {e:arealdensity/totalarealdensity for e,arealdensity in arealdensity.items()}
        
        layer1 = compoundfromlist.CompoundFromList(massfractions.keys(),massfractions.values(),types.fraction.mass,w.density,name=name)
        
        material = layer1
        thickness = windowthickness*1e-7
    else:
        elements = [compoundfromlist.CompoundFromList([e],[1],types.fraction.mole,0,name=e) for e in elements]
        layer1 = mixture.Mixture(elements,arealdensity,types.fraction.mass,name=name)

        layer2 = compoundfromname.compoundfromname("silicon nitride")

        material = [layer1,layer2]
        thickness = [filmthickness*1e-7,windowthickness*1e-7]
        
    return material,thickness,attenuators


class AXOID21_1(Standard):
    aliases = ["RF7-200-S2371-03"]
    
    def __init__(self,**kwargs):
        name = "RF7-200-S2371-03"
        elements = ["Pb","La","Pd","Mo","Cu","Fe","Ca"]
        arealdensity = [7.7,9,1.9,0.9,2.4,4,11.4] # ng/mm^2
        windowthickness = 200 # nm
        filmthickness = kwargs.pop("filmthickness",None) # nm
        material,thickness,attenuators = axo(name,elements,arealdensity,windowthickness,filmthickness)

        for k in attenuators:
            kwargs["geometry"].addattenuator(*k)

        super(AXOID21_1,self).__init__(material=material,thickness=thickness,**kwargs)


class AXOID21_2(Standard):
    aliases = ["RF8-200-S2454-17"]
    
    def __init__(self,**kwargs):
        name = "RF8-200-S2454-17"
        elements = ["Pb","La","Pd","Mo","Cu","Fe","Ca"]
        arealdensity = [6.3,7.6,2.3,0.7,2.6,4.1,25.1] # ng/mm^2
        windowthickness = 200 # nm
        filmthickness = kwargs.pop("filmthickness",None) # nm
        material,thickness,attenuators = axo(name,elements,arealdensity,windowthickness,filmthickness)
        
        for k in attenuators:
            kwargs["geometry"].addattenuator(*k)
            
        super(AXOID21_2,self).__init__(material=material,thickness=thickness,**kwargs)


class AXOID16b_1(Standard):
    aliases = ["RF8-200-S2453"]
    
    def __init__(self,**kwargs):
        name = "RF8-200-S2453"
        elements = ["Pb","La","Pd","Mo","Cu","Fe","Ca"]
        arealdensity = [5.9,10.3,1.2,.7,2.4,3.9,20.3] # ng/mm^2
        windowthickness = 200 # nm
        filmthickness = kwargs.pop("filmthickness",None) # nm
        material,thickness,attenuators = axo(name,elements,arealdensity,windowthickness,filmthickness)
        
        for k in attenuators:
            kwargs["geometry"].addattenuator(*k)
            
        super(AXOID16b_1,self).__init__(material=material,thickness=thickness,**kwargs)

factory = Standard.factory

