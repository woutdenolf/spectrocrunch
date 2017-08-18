# -*- coding: utf-8 -*-
#
#   Copyright (C) 2016 European Synchrotron Radiation Facility, Grenoble, France
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

from __future__ import division
import pyparsing
import math
import operator
import numpy as np
import scipy.ndimage.filters as filters
import h5py

import ..io.nexus as nexus

class StringParser(object):
    '''
    Most of this code comes from the fourFn.py pyparsing example
    
    http://stackoverflow.com/questions/2371436/evaluating-a-mathematical-expression-in-a-string
    '''
    
    def pushFirst(self, strg, loc, toks ):
        self.exprStack.append( toks[0] )
        
    def pushUMinus(self, strg, loc, toks ):
        if toks and toks[0]=='-': 
            self.exprStack.append( 'unary -' )
            
    def __init__(self,dtype=np.float32):
        """
        expop   :: '^'
        multop  :: '*' | '/'
        addop   :: '+' | '-'
        integer :: ['+' | '-'] '0'..'9'+
        atom    :: PI | E | real | fn '(' expr ')' | '(' expr ')'
        factor  :: atom [ expop factor ]*
        term    :: factor [ multop factor ]*
        expr    :: term [ addop term ]*
        """
        self.exprStack = []
        self.dtype = dtype

        point = pyparsing.Literal( "." )
        
        e     = pyparsing.CaselessLiteral( "E" )
        pi    = pyparsing.CaselessLiteral( "PI" )
        
        fnumber = pyparsing.Combine( pyparsing.Word( "+-"+pyparsing.nums, pyparsing.nums ) + 
                           pyparsing.Optional( point + pyparsing.Optional( pyparsing.Word( pyparsing.nums ) ) ) +
                           pyparsing.Optional( e + pyparsing.Word( "+-"+pyparsing.nums, pyparsing.nums ) ) )
                           
        ident = pyparsing.Word(pyparsing.alphas, pyparsing.alphas+pyparsing.nums+"_$")    
        var = pyparsing.Literal( "var_" ).suppress() +  pyparsing.Word(pyparsing.alphas)     

        plus  = pyparsing.Literal( "+" )
        minus = pyparsing.Literal( "-" )
        mult  = pyparsing.Literal( "*" )
        div   = pyparsing.Literal( "/" )
        lpar  = pyparsing.Literal( "(" ).suppress()
        rpar  = pyparsing.Literal( ")" ).suppress()
        addop  = plus | minus
        multop = mult | div
        expop = pyparsing.Literal( "^" )
        
        expr = pyparsing.Forward()
        atom = ((pyparsing.Optional(pyparsing.oneOf("- +")) +
                 (pi|e|fnumber|ident+lpar+expr+rpar|var).setParseAction(self.pushFirst))
                | pyparsing.Optional(pyparsing.oneOf("- +")) + pyparsing.Group(lpar+expr+rpar)
                ).setParseAction(self.pushUMinus)       

        # by defining exponentiation as "atom [ ^ factor ]..." instead of 
        # "atom [ ^ atom ]...", we get right-to-left exponents, instead of left-to-right
        # that is, 2^3^2 = 2^(3^2), not (2^3)^2.
        factor = pyparsing.Forward()
        factor << atom + pyparsing.ZeroOrMore( ( expop + factor ).setParseAction( self.pushFirst ) )

        term = factor + pyparsing.ZeroOrMore( ( multop + factor ).setParseAction( self.pushFirst ) )
        expr << term + pyparsing.ZeroOrMore( ( addop + term ).setParseAction( self.pushFirst ) )
        # addop_term = ( addop + term ).setParseAction( self.pushFirst )
        # general_term = term + ZeroOrMore( addop_term ) | OneOrMore( addop_term)
        # expr <<  general_term       
        self.bnf = expr
        
        # map operator symbols to corresponding arithmetic operations
        epsilon = 1e-12
        
        self.opn = { "+" : operator.add,
                "-" : operator.sub,
                "*" : operator.mul,
                "/" : operator.truediv,
                "^" : operator.pow }
                
        self.fn  = { "sin" : math.sin,
                "cos" : math.cos,
                "tan" : math.tan,
                "abs" : abs,
                "trunc" : lambda a: int(a),
                "round" : round,
                "sgn" : lambda a: abs(a)>epsilon and cmp(a,0) or 0,
                "max": np.nanmax,
                "min": np.nanmin,
                "mean": np.nanmean,
                "median": np.nanmedian if hasattr(np, 'nanmedian') else np.median,
                "ln": np.log,
                "sobel": lambda a: filters.sobel(a, axis=-1, output=None, mode='constant', cval=np.nan)}

        self.variables = {}

    def evaluateStack(self, s ):
        op = s.pop()
        if op == 'unary -':
            return -self.evaluateStack( s )
        if op in "+-*/^":
            op2 = self.evaluateStack( s )
            op1 = self.evaluateStack( s )
            return self.opn[op]( op1, op2 )
        elif op == "PI":
            return math.pi # 3.1415926535
        elif op == "E":
            return math.e  # 2.718281828
        elif op in self.fn:
            return self.fn[op]( self.evaluateStack( s ) )
        elif op[0].isalpha():
            if op in self.variables:
                var = self.variables[op]
                if isinstance(var,h5py.Group):
                    return var[var.attrs["signal"]][:]
                elif isinstance(var,h5py.Dataset):
                    return var[:]
                else:
                    return var
            return 0
        else:
            return self.dtype( op )
            
    def eval(self,num_string,parseAll=True,variables={}):
        self.exprStack=[]
        self.variables = variables
        results=self.bnf.parseString(num_string,parseAll)
        
        val=self.evaluateStack( self.exprStack[:] )
        return val

def evaluate_entire(operation,fin,varargs,fixedargs,retstacks):
    """ Evaluate expression on variable and fixed arguments
    Args:
        operation(dict)
        fin(nexus.File)
        varargs(dict)
        fixedargs(dict)
        retstacks(list(h5py.Group))
    """
    expression = operation["value"]

    mathparser = StringParser()
    ncalc = len(varargs)

    # Read fixed variables
    fixedvariables = {}
    for k in fixedargs:
        fixedvariables[k] = fin[fixedargs[k]]
        if isinstance(fixedvariables[k],h5py.Group):
            fixedvariables[k] = fixedvariables[k][fixedvariables[k].attrs["signal"]][:]
        if isinstance(fixedvariables[k],h5py.Dataset):
            fixedvariables[k] = fixedvariables[k][:]

    for i in range(ncalc):
        # Variables in expression
        variables = {}
        for k in varargs[i]:
            variables[k] = fin[varargs[i][k]]
        for k in fixedargs:
            variables[k] = fixedvariables[k]

        # Evaluate expression
        data = mathparser.eval(expression,variables=variables)

        # Add data to the NXdata group
        dset = nexus.createNXdataSignal(retstacks[i],data=data,chunks = True)

def evaluate_sliced(operation,fin,varargs,fixedargs,retstacks):
    """ Evaluate expression on variable and fixed arguments, slice by slice
    Args:
        operation(dict)
        fin(nexus.File)
        varargs(dict)
        fixedargs(dict)
        retstacks(list(h5py.Group))
    """
    expression = operation["value"]
    stackdim = operation["stackdim"]

    # Loop over variable arguments
    mathparser = StringParser()
    ncalc = len(varargs)
    for i in range(ncalc):
        # Variables in expression
        variables = {}
        for k in varargs[i]:
            variables[k] = fin[varargs[i][k]]
        for k in fixedargs:
            variables[k] = fin[fixedargs[k]]

        # Allocate space  
        shape = [0,0,0]
        v = np.int8(1)
        for k in variables:
            if isinstance(variables[k],h5py.Group):
                variables[k] = variables[k][variables[k].attrs["signal"]]

            if hasattr(variables[k],'shape'):
                s = variables[k].shape
                for j in range(min(len(s),3)):
                    shape[j] = max(shape[j],s[j])
                v *= variables[k][0,0,0]
            else:
                v *= variables[k]

        dset = nexus.createNXdataSignal(retstacks[i],shape=shape,dtype=type(v),chunks = True)

        # Evaluate expression
        for j in range(shape[stackdim]):
            # New (sliced) variables
            if stackdim==0:
                sliced = {k:variables[k][j,...] for k in variables}
            elif stackdim==1:
                sliced = {k:variables[k][:,j,:] for k in variables}
            else:
                sliced = {k:variables[k][...,j] for k in variables}

            data = mathparser.eval(expression,variables=sliced)

            if stackdim==0:
                dset[j,...] = data
            elif stackdim==1:
                dset[:,j,:] = data
            else:
                dset[...,j] = data

def evaluate(operation,fin,varargs,fixedargs,retstacks):
    if "sliced" in operation:
        if operation["sliced"]:
            evaluate_sliced(operation,fin,varargs,fixedargs,retstacks)
            return
    evaluate_entire(operation,fin,varargs,fixedargs,retstacks)

    
