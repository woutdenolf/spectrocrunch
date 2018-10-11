# -*- coding: utf-8 -*-
#
#   Copyright (C) 2016 European Synchrotron Radiation Facility, Grenoble, France
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


from ..patch.xraylib import xraylib
from ..patch.pint import ureg
from ..utils.hashable import Hashable
from . import element

import numpy as np

class Interaction(Hashable):

    def __init__(self,name,energy,prob):
        self._name = name
        self._energy = energy
        self._prob = prob

    @property
    def energy(self):
        return self._energy
        
    def _cmpkey(self):
        """For comparing
        """
        return str(self)
        
    def _sortkey(self):
        """For sorting
        """
        return self.energy
        
    def _stringrepr(self):
        """Unique representation of an instance
        """
        return self._name

class InteractionSource(Interaction):

    def __init__(self,energy,index):
        name = "Source-{}".format(index)
        prob = 1
        super(InteractionSource,self).__init__(name,energy,prob)

class InteractionFluo(Interaction):

    def __init__(self,el,shell,line):
        """
        Args:
            el(element or num or str):
            shell(Shell):
            line(FluoLine):
        """
        if not isinstance(el,element.Element):
            el = element.Element(el)
        name = "{}-{}".format(el,line)
        energy = line.energy(el.Z)
        prob = shell.fluoyield(el.Z)*line.radrate(el.Z)
        
        super(InteractionFluo,self).__init__(name,energy,prob)

class InteractionElScat(Interaction):

    def __init__(self,source):
        name = "RScat({})".format(source)
        prob = 1
        super(InteractionElScat,self).__init__(name,source.energy,prob)
        
class InteractionInelScat(Interaction):
        
    def __init__(self,source,theta):
        name = "CScat({})".format(source)
        prob = 1
        self.theta = theta # scattering angle
        
        super(InteractionInelScat,self).__init__(name,source.energy,prob)

    @property
    def energy(self):
        if self.theta==0:
            return self.energy
        delta = ureg.Quantity(1-np.cos(np.radians(self.theta)),"1/(m_e*c^2)").to("1/keV","spectroscopy").magnitude
        return self._energy/(1+self._energy*delta)
        

