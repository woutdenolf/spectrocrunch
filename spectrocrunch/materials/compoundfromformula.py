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

from . import compound
from .types import fractionType
import pyparsing as pp

class FormulaParser(object):

    def __init__(self):
        lpar  = pp.Literal( "(" ).suppress()
        rpar  = pp.Literal( ")" ).suppress()

        element = pp.Combine(pp.Word(pp.srange("[A-Z]"),exact=1)+pp.Optional(pp.Word(pp.srange("[a-z]"),max=1)))
        integer = pp.Word(pp.nums)
        point = pp.Literal(".")
        fnumber = pp.Combine(integer+pp.Optional(point+pp.Optional(integer))) | pp.Combine(point+integer)

        self.formula = pp.Forward()
        atom = element | pp.Group(lpar+self.formula+rpar)
        self.formula << pp.OneOrMore(pp.Group(atom+pp.Optional(fnumber, default="1")))
        self.elements = {}

    def parseresult(self,result,mult):
        """
        Args:
            result(pp.ParseResults): pyparsing result
            result(str): element or multiplier
            mult(float): multiplier
        """
        if isinstance(result,pp.ParseResults):
            if isinstance(result[-1],str):
                if not result[-1].isalpha():
                    mult *= float(result[-1])
            for r in result:
                self.parseresult(r,mult)
        elif result[-1].isalpha():
            if result in self.elements:
                self.elements[result] += mult
            else:
                self.elements[result] = mult
        else:
            pass

    def eval(self,formula):
        self.elements = {}
        result = self.formula.parseString(formula)
        self.parseresult(result,1.)
        return self.elements.keys(),self.elements.values()

class CompoundFromFormula(compound.Compound):
    """Interface to a compound defined by a chemical formula
    """

    def __init__(self,formula,density=None,name=None):
        """
        Args:
            formula(str): chemical formula
            density(float): compound density (g/cm^3) 
            name(Optional[str]): compound name
        """

        p = FormulaParser()
        elements, mults = p.eval(formula)
        if name is None:
            name = formula
        super(CompoundFromFormula,self).__init__(elements,mults,fractionType.mole,density,name=name)


