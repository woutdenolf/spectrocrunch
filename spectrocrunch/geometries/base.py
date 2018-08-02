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

import numpy as np
from ..common import units
from ..common import instance
from ..math import noisepropagation

class Base(object):

    def __init__(self,detector=None,source=None,atmosphere=None):
        self.detector = detector
        self.source = source
        self.atmosphere = atmosphere

    @property
    def detector(self):
        return self._detector
    
    @detector.setter
    def detector(self,value):
        self._detector = value
        if self._detector is not None:
            self._detector.geometry = self
    
    def __getattr__(self,attr):
        try:
            return getattr(self.detector,attr)
        except AttributeError:
            try:
                return getattr(self.source,attr)
            except AttributeError:
                raise AttributeError("'{}' object has no attribute '{}'".format(self.__class__.__name__,attr))
                
    def __str__(self):
        return "{}\n{}".format(self.source,self.detector)

    def addtofisx(self,setup,cfg):
        self.detector.addtofisx(setup,cfg)
    
    def addtopymca(self,setup,cfg):
        self.detector.addtopymca(setup,cfg)
        self.source.addtopymca(setup,cfg)
        
    def loadfrompymca(self,setup,cfg):
        self.detector.loadfrompymca(setup,cfg)
        self.source.loadfrompymca(setup,cfg)
        
        
class FlatSample(Base):

    def __init__(self,anglein=None,angleout=None,azimuth=0.,**kwargs):
        """
        Args:
            anglein(num): angle (deg) between primary beam and sample surface
            angleout(num): angle (deg) between detector and sample surface
            azimuth(num): angle (deg) between the source-detector plane and the polarization plane
        """
        self.anglein = anglein # deg
        self.angleout = angleout # deg
        self.azimuth = azimuth # deg
        
        super(FlatSample,self).__init__(**kwargs)

    @property
    def reflection(self):
        return self.angleout>0
        
    @property
    def cosnormin(self):
        # angle with surface normal (pointing inwards)
        return np.cos(np.radians(90-self.anglein))
    
    @property
    def cosnormout(self):
        # angle with surface normal (pointing inwards)
        return np.cos(np.radians(90+self.angleout))
    
    @property
    def scatteringangle(self):
        return self.anglein + self.angleout
    
    def xrayspectrumkwargs(self):
        return {"polar":np.radians(self.scatteringangle),"azimuth":np.radians(self.azimuth)}
        
    def __str__(self):
        return "{}\nGeometry:\n In = {} deg\n Out = {} deg ({})\n Azimuth = {} deg".format(\
                        super(FlatSample,self).__str__(),\
                        self.anglein,self.angleout,\
                        "reflection" if self.reflection else "transmission",\
                        self.azimuth)

    def addtofisx(self,setup,cfg):
        # When self.angleout<0: works only for a single layer
        setup.setGeometry(self.anglein, abs(self.angleout))
        super(FlatSample,self).addtofisx(setup,cfg)


class SolidAngle(FlatSample):

    def __init__(self,solidangle=None,**kwargs):
        self.solidangle = solidangle # srad
        super(SolidAngle,self).__init__(**kwargs)

    def __str__(self):
        if self.solidangle is None:
            return super(SolidAngle,self).__str__()
        else:
            return "{}\n Solid angle = 4*pi*{} srad".format(super(SolidAngle,self).__str__(),self.solidangle/(4*np.pi))
        
        
class Centric(FlatSample):

    def __init__(self,distance=None,**kwargs):
        """
        Args:
            distance(num): distance (cm) to target
        """
        self.distance = distance
        super(Centric,self).__init__(**kwargs)

    @property
    def distance_rv(self):
        return self._distance
    
    @property
    def distance(self):
        return noisepropagation.E(self.distance_rv)
        
    @distance.setter
    def distance(self,value):
        if value is None:
            self._distance = None
        else:
            self._distance = units.Quantity(distance,"cm")
        
    def calibrate_distance_manually(self):
        self.distance = distance
    
    @property
    def solidangle(self):
        return self.detector.solidangle_calc(activearea=self.detector.activearea,distance=self.distance)
        
    @solidangle.setter
    def solidangle(self,value):
        self.distance = self.detector.solidangle_calc(activearea=self.detector.activearea,solidangle=value)

    def __str__(self):
        if self.distance is None:
            return super(Centric,self).__str__()
        else:
            return "{}\n Distance = {:~}\n Solid angle = 4*pi*{} srad".format(super(Centric,self).__str__(),self.distance,self.solidangle/(4*np.pi))
        
    def addtopymca(self,setup,cfg): 
        super(Centric,self).addtopymca(setup,cfg)
        cfg["concentrations"]["distance"] = self.distance.to("cm").magnitude

    def loadfrompymca(self,setup,cfg): 
        super(Centric,self).loadfrompymca(setup,cfg)
        self.calibrate_distance_manually(cfg["concentrations"]["distance"])
         


