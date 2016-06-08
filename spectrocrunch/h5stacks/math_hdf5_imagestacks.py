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

from spectrocrunch.common.Enum import Enum
operationType = Enum(['expression','copy','crop','replace'])

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

def calccroproi(stack,nanval,stackdim):
    """Determine crop ROI so that no column or row consists of only nanval

    Args:
        stack(h5py.Group): NXdata group
        nanval(number)
        stackdim(int)

    Returns:
        tuple or None
    """

    dataset = stack[stack.attrs["signal"]]

    # Stack dimensions
    if stackdim==0:
        nimg,dim1,dim2 = dataset.shape
    elif stackdim==1:
        dim1,nimg,dim2 = dataset.shape
    else:
        dim1,dim2,nimg = dataset.shape

    # Mask (True = valid pixel)
    mask = np.ones((dim1,dim2),dtype=np.bool)
    for i in range(nimg):
        if stackdim==0:
            img = dataset[i,...]
        elif stackdim==1:
            img = dataset[:,i,:]
        else:
            img = dataset[...,i]
        if nanval is np.nan:
            mask &= np.isnan(img)==False
        else:
            mask &= img != nanval

    # Valid row and columns (not all False)
    indvalidrow = np.argwhere(mask.sum(axis=1))
    indvalidcol = np.argwhere(mask.sum(axis=0))
    arow = indvalidrow[0]
    brow = indvalidrow[-1]+1
    acol = indvalidcol[0]
    bcol = indvalidcol[-1]+1

    # Roi to keep
    if brow-arow == dim1 and bcol-acol == dim2:
        return None
    if stackdim==0:
        roi = ((0,nimg),\
               (arow,brow),\
               (acol,bcol))
    elif stackdim==1:
        roi = ((arow,brow),\
               (0,nimg),\
               (acol,bcol))
    else:
        roi = ((arow,brow),\
               (acol,bcol),\
               (0,nimg))
    return roi
    

def math_hdf5_imagestacks(filein,fileout,axes,operation,varargs,fixedargs,ret,extension="",overwrite=False,info=None,copygroups=None):
    """Perform some operations on hdf5 imagestacks

    Args:
        filein:
        fileout:
        axes(list(str)): list of NXdata groups
        operation(dict): type, value
        varargs(list(dict)|list(array-like)): 
        fixedargs(dict):
        ret:
        extension:
        overwrite:
        info:
        copygroups:
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

    # Create NXdata groups
    retstacks = [nexus.newNXdata(fout,s,extension) for s in ret]

    # Operation
    axesdata = []
    if operation["type"]==operationType.expression:
        expression = operation["value"]

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

            # Evaluate expression
            data = mathparser.eval(expression,variables=variables)

            # Add data to the NXdata group
            dset = nexus.createNXdataSignal(retstacks[i],data=data,chunks = True)

    elif operation["type"]==operationType.crop:
        reference = [s for s in varargs if s.endswith(operation["reference set"])]
        iref = varargs.index(reference[0])

        roi = calccroproi(fin[varargs[iref]],operation["nanval"],operation["stackdim"])

        # Add (cropped) data to the NXdata group
        if roi is None:
            for i in range(len(varargs)):
                grp = fin[varargs[i]]
                data = grp[grp.attrs["signal"]]
                dset = nexus.createNXdataSignal(retstacks[i],data=data,chunks = True)
        else:
            for i in range(len(varargs)):
                grp = fin[varargs[i]]
                data = grp[grp.attrs["signal"]][roi[0][0]:roi[0][1],roi[1][0]:roi[1][1],roi[2][0]:roi[2][1]]
                dset = nexus.createNXdataSignal(retstacks[i],data=data,chunks = True)

            axesdata = [None]*len(axes)
            for i in range(len(axes)):
                axesdata[i] = fin[axes[i]["fullname"]][roi[i][0]:roi[i][1]]
    elif operation["type"]==operationType.copy:
        for i in range(len(varargs)):
            grp = fin[varargs[i]]
            data = grp[grp.attrs["signal"]]
            dset = nexus.createNXdataSignal(retstacks[i],data=data,chunks = True)
    elif operation["type"]==operationType.replace:
        v1 = operation["v1"]
        v2 = operation["v2"]
        for i in range(len(varargs)):
            grp = fin[varargs[i]]
            data = grp[grp.attrs["signal"]].value
            if v1 is np.nan:
                data[np.isnan(data)] = v2
            else:
                data[data==v1] = v2
            dset = nexus.createNXdataSignal(retstacks[i],data=data,chunks = True)

    # New axes
    if len(axesdata)==0:
        # Copy of the old axes
        if bsame:
            axesdata = [fin[a["fullname"]] for a in axes]
        else:
            axesdata = [fin[a["fullname"]][:] for a in axes]
    retaxes = nexus.newaxes(fout,axes,axesdata,extension)

    # Link groups to axes
    nexus.linkaxes(fout,retaxes,retstacks)

    # Keep only names
    retstacks = [dset.name for dset in retstacks]

    # Return
    fin.close()
    if not bsame:
        fout.close()

    return retstacks, retaxes

def fluxnorm_hdf5_imagestacks(filein,fileout,axes,I0stack,stacks,retstack,overwrite=False,info=None,copygroups=None):
    operation = {"type":operationType.expression,"value":"var_a/var_b"}

    fixedargs = {"b":I0stack}

    nI = len(stacks)
    varargs = [None]*nI
    for i in range(nI):
        varargs[i] = {"a":stacks[i]}

    return math_hdf5_imagestacks(filein,fileout,axes,operation,varargs,fixedargs,retstack,extension="norm",overwrite=overwrite,info=info,copygroups=copygroups)

def copy_hdf5_imagestacks(filein,fileout,axes,stacks,retstack,overwrite=False,info=None,copygroups=None):
    operation = {"type":operationType.copy}
    return math_hdf5_imagestacks(filein,fileout,axes,operation,stacks,{},retstack,extension="copy",overwrite=overwrite,info=info,copygroups=copygroups)

def crop_hdf5_imagestacks(filein,fileout,axes,stacks,retstack,cropinfo,overwrite=False,info=None,copygroups=None):
    operation = cropinfo
    operation["type"] = operationType.crop
    return math_hdf5_imagestacks(filein,fileout,axes,operation,stacks,{},retstack,extension="crop",overwrite=overwrite,info=info,copygroups=copygroups)

def replacevalue_hdf5_imagestacks(filein,fileout,axes,stacks,retstack,orgvalue,newvalue,overwrite=False,info=None,copygroups=None):
    operation = {"type":operationType.replace,"v1":orgvalue,"v2":newvalue}
    return math_hdf5_imagestacks(filein,fileout,axes,operation,stacks,{},retstack,extension="replace",overwrite=overwrite,info=info,copygroups=copygroups)


