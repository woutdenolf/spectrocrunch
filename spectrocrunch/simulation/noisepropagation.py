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

from ..common.Enum import Enum

from uncertainties import unumpy
import numpy as np

class process(object):
    def expectedvalue(self):
        raise NotImplementedError()

    def variance(self):
        raise NotImplementedError()

class bernouilli(process):
    def __init__(self,probsuccess):
        self.probsuccess = probsuccess
        super(bernouilli,self).__init__()

    def expectedvalue(self):
        return self.probsuccess

    def variance(self):
        return self.probsuccess*(1-self.probsuccess)

class poisson(process):
    def __init__(self,gain):
        self.gain = gain
        super(poisson,self).__init__()

    def expectedvalue(self):
        return self.gain

    def variance(self):
        return self.gain

def propagate(N,process):
    """Sum of a random number (N) of independent random variables (X_i) with the same probability mass function
       
    Args:
        N(uncertainties.unumpy.uarray): incomming number of photons with uncertainties
        process(process): statistical process

    Returns:
        uncertainties.unumpy.uarray: X_1 + ... + X_N
    """
    EN = unumpy.nominal_values(N)
    VARN = unumpy.std_devs(N)**2

    EX = process.expectedvalue()
    VARX = process.variance()

    EY = EN*EX
    VARY = VARN*EX*EX + VARX*EN
        
    return unumpy.uarray(EY,np.sqrt(VARY))

