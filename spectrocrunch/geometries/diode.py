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

    def __init__(self,anglein=None,angleout=None):
        """
        Args:
            anglein(num): angle (deg) between primary beam and surface normal (pointing inwards)
            angleout(num): angle (deg) between fluorescene path to detector and surface normal (pointing inwards)
        """
        
        self.anglein = float(anglein) # deg
        self.angleout = float(angleout) # deg

    @property
    def cosanglein(self):
        return np.cos(np.radians(self.anglein))
    
    @property
    def cosangleout(self):
        return np.cos(np.radians(angleout))
    
    def __str__(self):
        return "In = {} deg\n Out = {} deg".format(self.distance,self.anglein,self.angleout)
        
factory = Geometry.factory

