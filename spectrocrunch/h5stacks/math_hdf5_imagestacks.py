# -*- coding: utf-8 -*-
#
#   Copyright (C) 2015 European Synchrotron Radiation Facility, Grenoble, France
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
import spectrocrunch.io.nexus as nexus

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
            
    def __init__(self):
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
                "median": np.nanmedian,
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
                else:
                    return var
            return 0
        else:
            return float( op )
            
    def eval(self,num_string,parseAll=True,variables={}):
        self.exprStack=[]
        self.variables = variables
        results=self.bnf.parseString(num_string,parseAll)
        
        val=self.evaluateStack( self.exprStack[:] )
        return val

def math_hdf5_imagestacks(filein,fileout,axes,expression,varargs,fixedargs,ret,extension="",overwrite=False,info=None,copygroups=None):
    """Perform some operations on hdf5 imagestacks

    Args:
        varargs(list(dist)): list of dictionaries
        fixedargs(dict): dictionary

    Returns:
        tuple(list,list)
    """

    bsame = filein==fileout

    # Open files
    if bsame:
        fin = nexus.File(filein,mode='r+')
        fout = fin
    
        if info is not None:
            nexus.addinfogroup(fout,"math",info)
    else:
        fin = nexus.File(filein,mode='r')
        fout = nexus.File(fileout,mode='w' if overwrite else 'a')
        extension = ""

        if info is not None:
            nexus.copyaddinfogroup(fin,fout,"math",info)

        if copygroups is not None:
            for grp in copygroups:
                fin.copy(fin[grp],fout)
                
    # New axes (just a copy, for clarity because other processes like alignment do change the axes)
    if bsame:
        axesdata = [fin[a["fullname"]] for a in axes]
    else:
        axesdata = [fin[a["fullname"]][:] for a in axes]
    retaxes = nexus.newaxes(fout,axes,axesdata,extension)

    # Create NXdata groups
    retstack = [nexus.newNXdata(fout,s,extension) for s in ret]

    # Loop over variable arguments
    mathparser = StringParser()
    ncalc = len(varargs)
    for icalc in range(ncalc):
        # Variables in expression
        variables = {}
        for k in varargs[icalc]:
            variables[k] = fin[varargs[icalc][k]]
        for k in fixedargs:
            variables[k] = fin[fixedargs[k]]

        # Evaluate expression
        data = mathparser.eval(expression,variables=variables)

        # Add data to the NXdata group
        dset = nexus.createNXdataSignal(retstack[icalc],data=data,chunks = True)

    # Link groups to axes
    nexus.linkaxes(fout,retaxes,retstack)

    # Keep only names
    retstack = [dset.name for dset in retstack]

    # Return
    fin.close()
    if not bsame:
        fout.close()

    return retstack, retaxes

def fluxnorm_hdf5_imagestacks(filein,fileout,axes,I0stack,stacks,retstack,overwrite=False,info=None,copygroups=None):
    fixedargs = {"b":I0stack}

    nI = len(stacks)
    varargs = [None]*nI
    for i in range(nI):
        varargs[i] = {"a":stacks[i]}

    return math_hdf5_imagestacks(filein,fileout,axes,"var_a/var_b",varargs,fixedargs,retstack,extension="norm",overwrite=overwrite,info=info,copygroups=copygroups)

def copy_hdf5_imagestacks(filein,fileout,axes,stacks,retstack,overwrite=False,info=None,copygroups=None):
    fixedargs = {}

    nI = len(stacks)
    varargs = [None]*nI
    for i in range(nI):
        varargs[i] = {"a":stacks[i]}

    return math_hdf5_imagestacks(filein,fileout,axes,"var_a",varargs,fixedargs,retstack,overwrite=overwrite,info=info,copygroups=copygroups)

