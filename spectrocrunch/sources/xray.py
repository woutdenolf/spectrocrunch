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


class Source(with_metaclass(object)):

    POLARIZATION = Enum(['none','linear','elliptical'])

    def __init__(self,polarization=None,Plinear=None,delta=None):
        """
        Args:
            polarization(): 
            Plinear(num): linear degree of polarization of the source
            delta(num): component phase retardation in deg (0 -> linear polarization, else elliptical) 
        """
        
        self.polarization = polarization
        
        if self.polarization == self.POLARIZATION.none:
            self.Plinear = 0
            self.delta = 90
        elif self.polarization == self.POLARIZATION.linear:
            self.Plinear = 1
            
        self.Plinear = float(Plinear)
        self.delta = float(delta)

    def __str__(self):
        return "Source:\n Linear degree of polarization = {} \n Phase retardation = {} deg".format(self.Plinear,self.delta)

    def addtopymca(self,setup,cfg):
        pass
    
    def loadfrompymca(self,setup,cfg):
        pass
            
    def K(self):
        if self.polarization == self.POLARIZATION.none:
            func = lambda theta,phi: (1+np.cos(theta)**2.)/2.
        elif self.polarization == self.POLARIZATION.linear:
            func = lambda theta,phi: (1+np.cos(theta)**2.-np.cos(2*phi)*np.sin(theta)**2.)/2.
        elif self.polarization == self.POLARIZATION.elliptical:
            func = lambda theta,phi: (1+np.cos(theta)**2.-(self.Plinear*np.cos(2*phi)-np.sqrt(1-self.Plinear**2.)*np.sin(2*phi)*self.cosdelta*np.sin(theta)**2.))/2.
            
        return func

    @property
    def cosbeta(self):
        return np.sqrt((1+self.Plinear)/2.)

    @property
    def sinbeta(self):
        return np.sqrt((1-self.Plinear)/2.)

    @property
    def sin2beta(self):
        return np.sqrt(1-self.Plinear**2.)/2.

    @property
    def cos2beta(self):
        return self.Plinear
    
    @property
    def tan2beta(self):
        return np.sqrt(1./self.Plinear**2.-1)/2.
            
    @property
    def cosdelta(self):
        return np.cos(np.radians(self.delta))

    @property
    def sindelta(self):
        return np.sin(np.radians(self.delta))
        
    def plotellipse(self,meanabsphasor):
        if self.polarization == POLARIZATION.none:
            raise RuntimeError("Unpolarized radiation has a random electric field vector")
    
        semimajor = meanabsphasor*np.sqrt((1+np.sqrt(1-(self.sin2beta*self.sindelta)**2))/2.)
        semiminor = meanabsphasor*np.sqrt((1-np.sqrt(1-(self.sin2beta*self.sindelta)**2))/2.)
        chi = np.arctan(self.tan2beta*self.cosdelta)/2.
        
        p = np.linspace(0, 2*np.pi, len(t))
        x = semimajor * np.cos(p) * np.cos(chi) - semiminor * np.sin(p) * np.sin(chi)
        y = semimajor * np.cos(p) * np.sin(chi) + semiminor * np.sin(p) * np.cos(chi)
        plt.plot(x,y)

        x = [0,semimajor*np.cos(chi)]
        y = [0,semimajor*np.sin(chi)]
        plt.plot(x,y)
        
        x = [0,semiminor*np.cos(chi+np.pi/2)]
        y = [0,semiminor*np.sin(chi+np.pi/2)]
        plt.plot(x,y)
        
        plt.gca().set_aspect('equal', 'datalim')
    
    def ellipseeq(self,meanabsphasor):
        """ A.x^2 + B.x.y + C.y^2 = 1
        """
        A = 1/(meanabsphasor*self.sindelta*self.cosbeta)**2
        C = 1/(meanabsphasor*self.sindelta*self.sinbeta)**2
        B = -2*self.cosdelta/(meanabsphasor**2*self.sindelta**2*self.sinbeta*self.cosbeta)
        return A,B,C
        
        
class Synchrotron(Source):

    def __init__(self,**kwargs):
        #super(Synchrotron,self).__init__(polarization=Source.POLARIZATION.linear)
        super(Synchrotron,self).__init__(polarization=Source.POLARIZATION.elliptical,Plinear=0.95,delta=90)


class Tube(Source):

    def __init__(self,**kwargs):
        super(Tube,self).__init__(polarization=Source.POLARIZATION.none)


factory = Source.factory

