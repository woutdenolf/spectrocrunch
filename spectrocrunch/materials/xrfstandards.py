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

import numpy as np

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


class ThinFilmStandard(Standard):
    
    def __init__(self,arealdensity,substrate,substratethickness,name=None,filmthickness=None,**kwargs):
        arealdensity_film = dict(zip(map(element.Element,arealdensity.keys()),arealdensity.values()))

        if filmthickness is None:
            wfrac_substrate = substrate.elemental_massfractions()
            
            elfilm = set(arealdensity_film.keys())
            elsubstrate = set(wfrac_substrate.keys())
            if elfilm & elsubstrate:
                raise RuntimeError("Film and substrate cannot contain the same elements")
                
            kwargs["material"] = substrate
            kwargs["thickness"] = substratethickness
            if "extra" not in kwargs:
                kwargs["extra"] = [] 
            kwargs["extra"].extend(elfilm)
            
            self._arealdensities = arealdensity_film
            m = substrate.density*substratethickness
            self._arealdensities.update(dict(zip(wfrac_substrate.keys(),np.array(wfrac_substrate.values())*m)))
        else:
            ad = np.array(arealdensity_film.values())
            sad = ad.sum()
            density = sad/filmthickness
            massfractions = ad/sad
            film = compoundfromlist.CompoundFromList(arealdensity_film.keys(),massfractions,\
                                                     types.fraction.mass,density=density,name=name)
            kwargs["material"] = [film,substrate]
            kwargs["thickness"] = [filmthickness,substratethickness]
        
            self._arealdensities = None
        
        super(ThinFilmStandard,self).__init__(**kwargs)
        
    def arealdensity(self):
        if self._arealdensities is None:
            return super(ThinFilmStandard,self).arealdensity()
        else:
            return self._arealdensities


class AXOID21_1(ThinFilmStandard):
    aliases = ["RF7-200-S2371-03","id21_room"]
    
    def __init__(self,**kwargs):
        name = "RF7-200-S2371-03"
        
        # ng/mm^2 -> g/cm^2
        arealdensity = {"Pb": 7.7e-7,\
                        "La": 9.0e-7,\
                        "Pd": 1.9e-7,\
                        "Mo": 0.9e-7,\
                        "Cu": 2.4e-7,\
                        "Fe": 4.0e-7,\
                        "Ca":11.4e-7}
        filmthickness = kwargs.pop("filmthickness",None) # cm
        substrate = compoundfromname.compoundfromname("silicon nitride")
        substratethickness = 200e-7 # cm
        
        ultralene = compoundfromname.compoundfromname("ultralene")
        attenuators = [["SampleCover",ultralene,4e-4],\
                       ["BeamFilter0",ultralene,4e-4]]
        for k in attenuators:
            kwargs["geometry"].addattenuator(*k)

        super(AXOID21_1,self).__init__(arealdensity,substrate,substratethickness,name=name,filmthickness=filmthickness,**kwargs)


class AXOID21_2(ThinFilmStandard):
    aliases = ["RF8-200-S2454-17","id21_cryo"]
    
    def __init__(self,**kwargs):
        name = "RF8-200-S2454-17"

        # ng/mm^2 -> g/cm^2
        arealdensity = {"Pb": 6.3e-7,\
                        "La": 7.6e-7,\
                        "Pd": 2.3e-7,\
                        "Mo": 0.7e-7,\
                        "Cu": 2.6e-7,\
                        "Fe": 4.1e-7,\
                        "Ca":25.1e-7}
        filmthickness = kwargs.pop("filmthickness",None) # cm
        substrate = compoundfromname.compoundfromname("silicon nitride")
        substratethickness = 200e-7 # cm
        
        ultralene = compoundfromname.compoundfromname("ultralene")
        attenuators = [["SampleCover",ultralene,4e-4],\
                       ["BeamFilter0",ultralene,4e-4]]
        for k in attenuators:
            kwargs["geometry"].addattenuator(*k)

        super(AXOID21_2,self).__init__(arealdensity,substrate,substratethickness,name=name,filmthickness=filmthickness,**kwargs)
        
        
class AXOID16b_1(ThinFilmStandard):
    aliases = ["RF8-200-S2453"]
    
    def __init__(self,**kwargs):
        name = "RF8-200-S2453"

        # ng/mm^2 -> g/cm^2
        arealdensity = {"Pb": 5.9e-7,\
                        "La": 10.3e-7,\
                        "Pd": 1.2e-7,\
                        "Mo": 0.7e-7,\
                        "Cu": 2.4e-7,\
                        "Fe": 3.9e-7,\
                        "Ca":20.3e-7}
        filmthickness = kwargs.pop("filmthickness",None) # cm
        substrate = compoundfromname.compoundfromname("silicon nitride")
        substratethickness = 200e-7 # cm

        super(AXOID16b_1,self).__init__(arealdensity,substrate,substratethickness,name=name,filmthickness=filmthickness,**kwargs)
        

factory = Standard.factory

