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
from ..common.instance import isarray

import numpy as np

from uncertainties.core import Variable as RandomVariable

from uncertainties import unumpy

# Uncertainties: linear approximation to error propagation
#
#   Variable.std_dev():
#       sqrt(Da^2 + Db^2 + ...)
#
#   Variable.error_components():
#       {a:Da,b:Db,...}
#       Da = abs(da * sa)
#
#   Variable.derivatives:
#       {a:da,b:db,...}
#   
#   covariance_matrix([x,y,...]):
#       C[x,y] = sum( x.derivatives[a] * y.derivatives[a] * sa^2 + 
#                     x.derivatives[b] * y.derivatives[b] * sb^2 + ... )
#

class Bernouilli(RandomVariable):
    def __init__(self,probsuccess,**kwargs):
        super(Bernouilli,self).__init__(probsuccess,np.sqrt(probsuccess*(1-probsuccess)),**kwargs)

class Poisson(RandomVariable):
    def __init__(self,gain,**kwargs):
        super(Poisson,self).__init__(gain,np.sqrt(gain),**kwargs)

def bernouilli(p):
    if isarray(p):
        return np.vectorize(Bernouilli,otypes=[object])(p)
    else:
        return Bernouilli(p)

def poisson(p):
    if isarray(p):
        return np.vectorize(Poisson,otypes=[object])(p)
    else:
        return Poisson(p)
        
def randomvariable(X,SX):
    if isarray(X):
        return np.vectorize(lambda x, s: RandomVariable(x, s), otypes=[object])(X,SX)
    else:
        return RandomVariable(X,SX)
        
def E(X):
    if isinstance(X,RandomVariable):
        return X.nominal_value
    else:
        return unumpy.nominal_values(X)
    
def S(X):
    if isinstance(X,RandomVariable):
        return X.std_dev
    else:
        return unumpy.std_devs(X)

def VAR(X):
    return S(X)**2
        
def SNR(X):
    return E(X)/S(X)

def NSR(X):
    return S(X)/E(X)
    
def RVAR(X):
    return VAR(X)/E(X)**2
    
def repeat(N,X):
    """Sum of a fixed number (N) of independent random variables (Xi) with the same probability mass function.

        ISSUE: loose correlation with X

    Args:
        N(num): number of repeats
        X(num|array): random variable
        
    """
    # Linear error propagation of X+X:
    #    VAR[X+X] = VAR[X] + VAR[X] + 2.COV[X,X]
    #             = VAR[X] + VAR[X] + 2.VAR[X]
    #             = 4.VAR[X]
    #   or equivalently
    #    VAR[2.X] = 2^2.VAR[X]
    #   or using the uncertainties package
    #    Y = N*X
    #
    # => SNR[n.X] = SNR[X]
    #
    # However we want repeats in the sense COV[X,X] = 0
    #    VAR[X+X] = VAR[X] + VAR[X]
    #   or equivalently
    #     propagate(ufloat(N,0),X)
    #   or using the uncertainties package
    #     see below
    #
    # => SNR[n.X] = SNR[X] * sqrt(n)
    
    return randomvariable(N*E(X),np.sqrt(N)*S(X))

def compound(N,X):
    """Sum of a random number (N) of independent random variables (X_i) with E[X_i]=E[X_j], VAR[X_i]=VAR[X_j]
       (which is weaker than saying the have the same pmf). This is called a "compound random variable".
       
       S. K. Ross, Introduction to Probability Models, Eleventh Edition, Academic Press, 2014, example 3.10 and 3.19.
       
       http://www.math.unl.edu/~sdunbar1/ProbabilityTheory/Lessons/Conditionals/RandomSums/randsum.shtml
       https://mathmodelsblog.wordpress.com/2010/01/17/an-introduction-to-compound-distributions/
       
        ISSUE: loose correlation with N and X
        
    Args:
        N(num|array): incomming number of photons with uncertainties
        X(num|array): probability mass function of X

    Returns:
        array: Y = X_1 + ... + X_N  (nX x nN)
    """

    # 1. pmf(X) = Bernouilli (e.g. X-ray transmission) ->  Binomial selection
    #       -> N a fixed number: pmf(Y) = Binomial
    #       -> pmf(N) = Poisson: pmf(Y) = Poisson (binomial selection theorem)
    # 2. pmf(X) = Poisson  (e.g. X-ray to VIS)  
    
    bexpand = isarray(X) or isarray(N)
    
    EN = E(N)
    VARN = VAR(N)
    EX = E(X)
    VARX = VAR(X)
    
    if bexpand:
        nN = np.asarray(N).shape
        if len(nN)==2:
            nN = nN[1]
        else:
            nN = int(np.product(nN))
        nX = np.asarray(X).size

        EN = np.broadcast_to(EN,[nX,nN])
        VARN = np.broadcast_to(VARN,[nX,nN])

        EX = np.broadcast_to(EX,[nN,nX]).T
        VARX = np.broadcast_to(VARX,[nN,nX]).T
    
    EY = EN*EX
    VARY = VARN*EX*EX + VARX*EN

    return randomvariable(EY,np.sqrt(VARY))


