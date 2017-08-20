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

import numpy as np
from uncertainties import unumpy
from uncertainties import ufloat

# bernouilli and poisson could be derived from uncertainties.Variable
# but I want to prevent the sqrt(VAR)**2 round-off errors.

class bernouilli(object):
    def __init__(self,probsuccess):
        self.probsuccess = probsuccess

    def expectedvalue(self):
        return self.probsuccess

    def variance(self):
        return self.probsuccess*(1-self.probsuccess)

class poisson(object):
    def __init__(self,gain):
        self.gain = gain

    def expectedvalue(self):
        return self.gain

    def variance(self):
        return self.gain

def randomvariable(X,SX):
    if hasattr(X,"__iter__"):
        return unumpy.uarray(X,SX)
    else:
        return ufloat(X,SX)

def repeat(N,X):
    """Sum of a fixed number (N) of independent random variables (Xi) with the same probability mass function.

    Args:
        N(num): number of repeats
        X(uncertainties.unumpy.uarray): random variable
        
    """
    # Normal error propagation assumes COV[Xi,Xj] = VAR[X]
    #   For example:
    #    VAR[X+X] = VAR[X] + VAR[X] + 2.COV[X,X]
    #             = VAR[X] + VAR[X] + 2.VAR[X]
    #             = 4.VAR[X]
    #   or equivalently
    #    VAR[2.X] = 2^2.VAR[X]
    #   or using the uncertainties package
    #    Y = N*X
    #
    # We will assume that COV[Xi,Xj] = 0
    #   For example:
    #    VAR[X+X] = VAR[X] + VAR[X]
    #   or equivalently
    #     propagate(N,X) with VARN=0
    #   or using the uncertainties package
    #     sum([copy.copy(X) for i in range(N)])

    if hasattr(X,"__iter__"):
        return unumpy.uarray(N*unumpy.nominal_values(X),np.sqrt(N)*unumpy.std_devs(X))
    else:
        return ufloat(N*X.n,np.sqrt(N)*X.s)
    
def propagate(N,X):
    """Sum of a random number (N) of independent random variables (X_i) with the same probability mass function.
       
       http://www.math.unl.edu/~sdunbar1/ProbabilityTheory/Lessons/Conditionals/RandomSums/randsum.shtml
       
    Args:
        N(uncertainties.unumpy.uarray): incomming number of photons with uncertainties
        X(pmf): probability mass function of X

    Returns:
        uncertainties.unumpy.uarray: Y = X_1 + ... + X_N
    """

    # 1. pmf(X) = Bernouilli (e.g. X-ray transmission) ->  Binomial selection
    #       -> N a fixed number: pmf(Y) = Binomial
    #       -> pmf(N) = Poisson: pmf(Y) = Poisson (binomial selection theorem)
    # 2. pmf(X) = Poisson  (e.g. X-ray to VIS)  

    EN = unumpy.nominal_values(N)
    VARN = unumpy.std_devs(N)**2

    EX = X.expectedvalue()
    VARX = X.variance()

    EY = EN*EX
    VARY = VARN*EX*EX + VARX*EN
    
    return randomvariable(EY,np.sqrt(VARY))

