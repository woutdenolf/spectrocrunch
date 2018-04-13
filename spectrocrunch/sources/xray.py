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

from ..common.classfactory import with_metaclass
from ..common.Enum import Enum
from . import polarization

import numpy as np
import matplotlib.pyplot as plt

class Source(with_metaclass(object)):

    def __init__(self,stokes=None):
        self.stokes = stokes
    
    def __str__(self):
        s = str(self.stokes).replace('\n','\n ')
        return "Source:\n {}".format(s)
    
    @property
    def thomson_K(self):
        return self.stokes.thomson_K
    
    def compton_K(self,energy):
        return self.stokes.compton_K(energy)
        
    def addtopymca(self,setup,cfg):
        pass
    
    def loadfrompymca(self,setup,cfg):
        pass
    
    
class Synchrotron(Source):

    def __init__(self,**polparams):
        if "intensity" not in polparams:
            polparams["intensity"] = 1 # W/m^2
        if "dop" not in polparams:
            polparams["dop"] = 0.7# degree of polarization (in [0,1])
        if "dolp" not in polparams:
            polparams["dolp"] = 0.9*polparams["dop"]# degree of linear polarization (in [0,dop])
        if "polangle" not in polparams:
            polparams["polangle"] = 0 # angle of polarization ellipse with respect to the horizontal direction (in [-90,90])
        if "handedness" not in polparams:
            polparams["handedness"] = "left" # above/below the plane
        stokes = polarization.Stokes.from_params(**polparams)
        super(Synchrotron,self).__init__(stokes=stokes)


class Tube(Source):

    def __init__(self,intensity=1):
        stokes = polarization.Stokes.from_params(intensity=intensity,dop=0,dolp=0,polangle=0,handedness="left")
        super(Tube,self).__init__(stokes=stokes)


factory = Source.factory

