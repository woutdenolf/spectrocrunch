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

class Base(object):

    def __init__(self,detector=None,source=None):
        self.detector = detector
        self.source = source
        
        self.detector.geometry = self
        self.source.geometry = self
    
    def __getattr__(self,attr):
        try:
            return getattr(self.detector,attr)
        except AttributeError:
            try:
                return getattr(self.source,attr)
            except AttributeError:
                raise AttributeError("'{}' object has no attribute '{}'".format(self.__class__.__name__,subject))
                
    def __str__(self):
        return "Source:\n{}\nDetector:\n{}".format(self.source,self.detector)

    def addtofisx(self,setup,cfg):
        self.detector.addtofisx(setup,cfg)
        
        
class FlatSample(Base):

    def __init__(self,anglein=None,angleout=None,**kwargs):
        """
        Args:
            anglein(num): angle (deg) between primary beam and sample surface
            angleout(num): angle (deg) between detector and sample surface
        """
        self.anglein = float(anglein) # deg
        self.angleout = float(angleout) # deg
        
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
    
    def __str__(self):
        return "{}\nGeometry:\n In = {} deg\n Out = {} deg ({})".format(\
                        super(FlatSample,self).__str__(),\
                        self.anglein,self.angleout,"reflection" if self.reflection else "transmission")

    def addtofisx(self,setup,cfg):
        # When self.angleout<0: works only for a single layer
        setup.setGeometry(self.anglein, abs(self.angleout))
        super(FlatSample,self).addtofisx(setup,cfg)


class Point(FlatSample):

    def __init__(self,azimuth=None,**kwargs):
        """
        Args:
            azimuth(num): angle (deg) between the source-detector plane and the polarization plane
        """
        
        self.azimuth = float(azimuth) # deg
        
        super(Point,self).__init__(**kwargs)

    def xrayspectrumkwargs(self):
        return {"polar":self.scatteringangle,"azimuth":self.azimuth}
        
        
