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
from . import base

import numpy as np

class Geometry(with_metaclass(base.Point)):

    def __init__(self,detectorposition=None,distancefunc=None,distanceifunc=None,\
                **kwargs):
        """
        Args:
            detectorposition(num): motor position in motor units
            distancefunc(callable): convert detectorposition to distance in cm
            distanceifunc(callable): inverse of distancefunc
        """
        
        self.detectorposition = float(detectorposition)
        if distancefunc is None or distanceifunc is None:
            distancefunc = lambda x:x
            distanceifunc = lambda x:x
        self.distancefunc = distancefunc
        self.distanceifunc = distanceifunc

        super(Geometry,self).__init__(**kwargs)

    @property
    def distance(self):
        """Sampel detector distance in cm
        """
        return self.distancefunc(self.detectorposition)
        
    @distance.setter
    def distance(self,value):
        self.detectorposition = self.distanceifunc(value)

    def __str__(self):
        return "{}\n Distance = {} cm".format(super(Geometry,self).__str__(),self.distance)
        
        
class sdd120(Geometry):

    def __init__(self,**kwargs):
        # blc10516 (April 2017)
        # detector position in mm
        distancefunc = lambda x: (x+60.5)/10
        distanceifunc = lambda x: x*10-60.5
        super(sdd120,self).__init__(anglein=62,angleout=49,azimuth=0,\
                        distancefunc=distancefunc,distanceifunc=distanceifunc,\
                        **kwargs)

class sdd90(Geometry):

    def __init__(self,**kwargs):
        # blc10516 (April 2017)
        # detector position in mm
        distancefunc = lambda x: (x+85.5)/10
        distanceifunc = lambda x: x*10-85.5
        super(sdd90,self).__init__(anglein=62,angleout=28,azimuth=0,\
                        distancefunc=distancefunc,distanceifunc=distanceifunc,\
                        **kwargs)

factory = Geometry.factory

