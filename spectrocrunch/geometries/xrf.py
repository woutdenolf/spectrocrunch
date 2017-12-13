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
import numpy as np

class Geometry(with_metaclass(object)):

    def __init__(self,anglein=None,angleout=None,detectorposition=None,distanceoffset=0):
        """
        Args:
            anglein(num): angle (deg) between primary beam and surface
            angleout(num): angle (deg) between fluorescene path to detector and surface
            detectorposition(num): motor position in cm
            distanceoffset(num): position offset in cm to get the true distance
        """
        
        self.anglein = float(anglein) # deg
        self.angleout = float(angleout) # deg
        self.detectorposition = float(detectorposition)
        self.distanceoffset = float(distanceoffset)
        self.reflection = self.angleout>0

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
    
    @property
    def distance(self):
        return self.detectorposition + self.distanceoffset
        
    @distance.setter
    def distance(self,value):
        self.detectorposition = value - self.distanceoffset

    def __str__(self):
        return "Distance = {} mm\n In = {} deg\n Out = {} deg ({})".format(self.distance,self.anglein,self.angleout,"reflection" if self.reflection else "transmission")
        
class sdd120(Geometry):

    def __init__(self,**kwargs):
        # blc10516 (April 2017)
        super(sdd120,self).__init__(anglein=62,angleout=49,distanceoffset=6.05,**kwargs)

class sdd90(Geometry):

    def __init__(self,**kwargs):
        # blc10516 (April 2017)
        super(sdd90,self).__init__(anglein=62,angleout=28,distanceoffset=8.55,**kwargs)

factory = Geometry.factory

