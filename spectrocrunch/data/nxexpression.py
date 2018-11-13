# -*- coding: utf-8 -*-
#
#   Copyright (C) 2018 European Synchrotron Radiation Facility, Grenoble, France
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
import numpy as np
import scipy.ndimage.filters as filters
import h5py
import re
import string
import collections
import logging

from . import nxregulargrid
from ..utils.integerbase import integerbase
from ..io import nxfs

logger = logging.getLogger(__name__)

class Task(nxregulargrid.Task):
    
    def _parameters_defaults(self):
        super(Task,self)._parameters_defaults()
        self._required_parameters('expression')
        parameters = self.parameters
        parameters['skip'] = parameters.get('skip',[])
        logger.info('Skip regex: {}'.format(parameters['skip']))
        logger.info('Expression: {}'.format(parameters['expression']))

    def _parameters_filter(self):
        return super(Task,self)._parameters_filter()+['expression','skip']
    
    @staticmethod
    def skipfunc(skipdict):
        method = skipdict.get('method','regex')
        if method=='regexparent':
            return lambda signal:re.match(skipdict['pattern'],signal.parent.name)
        elif method=='regex':
            return lambda signal:re.match(skipdict['pattern'],signal.name)
        else:
            return lambda signal:False
        
    def _prepare_process(self):
        super(Task,self)._prepare_process()
        self.mathparser = MathParser()
        self._parse_expression()
        self.skipfuncs = [self.skipfunc(skipdict) for skipdict in self.parameters['skip']]

    def _process_data(self,data):
        if self.variables:
            variables = {'a':data}
            for name,path in self.variables.items():
                with path.open() as dset:
                    variables[name] = dset[tuple(self.indexin)]
            data = self.mathparser.eval(self.expression,variables=variables)
        return data
            
    def _parse_expression(self):
        # For example: "{}*nanone({xmap_icr}/({arr_iodet}*{xmap_ocr})"
        #  expression = 'var_a*nanone(var_b/(var_c*var_d))'
        #  namemap = {'var_b':'xmap_icr','var_c':'arr_iodet','var_d':'xmap_ocr'}
        
        varnames = self._extract_variable_names()
        if '' not in varnames:
            raise ValueError('The expression does not contain {} which indicates the variable argument.')

        # Rename variable argument
        expression = self.parameters['expression']
        expression = expression.replace('{}','var_a')
        varnames.remove("")
        
        # Rename fixed arguments
        o = integerbase(digs = [c for c in string.ascii_lowercase])
        namemap = {}
        for i,oldname in enumerate(varnames,1):
            newname = o.int2base(i)
            expression = expression.replace('{{{}}}'.format(oldname),'var_{}'.format(newname))
            namemap[oldname] = newname
        self.expression = expression
        self.namemap = namemap
        
        # Map of variables
        varcount = collections.Counter()
        vardefault = {}
        for signal in self.grid.signals:
            varcount[signal.name] += 1
            vardefault[signal.name] = signal
        self.vardefault = {k:vardefault[k] for k,v in varcount.items() if v==1}

    def _prepare_signal(self,signal):
        if self._skip(signal):
            logger.info('saving (skip) {}'.format(signal.name))
            self.variables = {}
        else:
            logger.info('saving {}'.format(signal.name))
            # For example: 'xmap_icr':'var_c' -> 'var_c':.../detector00/xmap_icr
            self.variables = {}
            for oldname,newname in self.namemap.items():
                path = self.vardefault.get(oldname,None)
                if not path:
                    path = signal.parent[oldname]
                self.variables[newname] = path

    def _skip(self,signal):
        for func in self.skipfuncs:
            if func(signal):
                return True
        return False

    def _extract_variable_names(self,allowempty=True):
        expression = self.parameters['expression']
        if allowempty:
            return list(set(re.findall('\{(.*?)\}',expression)))
        else:
            return list(set(re.findall('\{([^\{\}]+)\}',expression)))
    
class MathParser(object):
    """
    Most of this code comes from the fourFn.py pyparsing example
    http://stackoverflow.com/questions/2371436/evaluating-a-mathematical-expression-in-a-string
    """
    
    def pushFirst(self, strg, loc, toks):
        self.exprStack.append(toks[0])
        
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
        
        self.opn = { "+" : np.add,
                "-" : np.subtract,
                "*" : np.multiply,
                "/" : np.true_divide,
                "^" : np.power }
                
        self.fn  = { "sin" : np.sin,
                "cos" : np.cos,
                "tan" : np.tan,
                "abs" : abs,
                "trunc" : lambda a: int(a),
                "round" : round,
                "sgn" : lambda a: abs(a)>epsilon and cmp(a,0) or 0, #TODO: or 0?
                "max": np.nanmax,
                "min": np.nanmin,
                "nanone": lambda a: self._fn_nan(a,1),
                "nanzero": lambda a: self._fn_nan(a,0),
                "mean": np.nanmean,
                "median": np.nanmedian if hasattr(np, 'nanmedian') else np.median,
                "ln": np.log,
                "sobel": lambda a: filters.sobel(a, axis=-1, output=None, mode='constant', cval=np.nan)}

        self.variables = {}

    @staticmethod
    def _fn_nan(a,v):
        a[np.isnan(a)] = v
        return a
        
    def evaluateStack(self,s):
        op = s.pop()
        if op == 'unary -':
            return -self.evaluateStack(s)
        if op in "+-*/^":
            op2 = self.evaluateStack(s)
            op1 = self.evaluateStack(s)
            return self.opn[op](op1,op2)
        elif op.upper() == "PI":
            return np.pi # 3.1415926535
        elif op.upper() == "E":
            return np.e  # 2.718281828
        elif op in self.fn:
            return self.fn[op](self.evaluateStack(s))
        elif op[0].isalpha():
            if op in self.variables:
                var = self.variables[op]
                if isinstance(var,nxfs.Path):
                    if var.isfile:
                        with open(var) as dset:
                            return dset[:]
                    else:
                        raise ValueError('{} must be a file'.format(var))
                elif isinstance(var,h5py.Dataset):
                    return var[:]
                else:
                    return var
            return 0
        else:
            return self.dtype(op)
            
    def eval(self,expression,parseAll=True,variables={}):
        self.exprStack = []
        self.variables = variables
        results = self.bnf.parseString(expression,parseAll)
        return self.evaluateStack(list(self.exprStack))
