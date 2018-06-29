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

from ..common import lut

class Optics(object):

    def __init__(self):
        self.lut = lut.LUT(default=1)
    
    def __str__(self):
        name = type(self).__name__
        s = '\n '.join("{} keV: {} %".format(k,v*100) for k,v in self.lut.table())
        if s:
            return "{}:\n transmission:\n {}".format(name,s)
        else:
            return "{}:\n transmission: 100%".format(name)
    
    def reset_transmission(self):
        if self.haslut():
            self.lut.clear(1)
      
    def transmission(self,energy):
        self.checklut()
        return self.lut(energy)
    
    def set_transmission(self,energy,transmission):
        self.checklut()
        self.lut.add(energy,transmission)

    def haslut(self):
        return hasattr(self,"lut")
    
    def checklut(self):
        if not self.haslut():
            raise RuntimeError("{} has no transmission lookup table.".format(type(self).__name__))
            
            
