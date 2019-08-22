# -*- coding: utf-8 -*-

from . import compound
from . import types
from ..utils import instance
import pyparsing as pp


class FormulaParser(object):

    def __init__(self):
        lpar = pp.Literal("(").suppress()
        rpar = pp.Literal(")").suppress()

        element = pp.Combine(pp.Word(
            pp.srange("[A-Z]"), exact=1)+pp.Optional(pp.Word(pp.srange("[a-z]"), max=1)))
        integer = pp.Word(pp.nums)
        point = pp.Literal(".")
        fnumber = pp.Combine(
            integer+pp.Optional(point+pp.Optional(integer))) | pp.Combine(point+integer)

        self.formula = pp.Forward()
        atom = element | pp.Group(lpar+self.formula+rpar)
        self.formula << pp.OneOrMore(
            pp.Group(atom+pp.Optional(fnumber, default="1")))
        self.elements = {}

    def parseresult(self, result, mult):
        """
        Args:
            result(pp.ParseResults): pyparsing result
            result(str): element or multiplier
            mult(float): multiplier
        """
        if isinstance(result, pp.ParseResults):
            if isinstance(result[-1], str):
                if not result[-1].isalpha():
                    mult *= float(result[-1])
            for r in result:
                self.parseresult(r, mult)
        elif result[-1].isalpha():
            if result in self.elements:
                self.elements[result] += mult
            else:
                self.elements[result] = mult
        else:
            pass

    def eval(self, formula):
        self.elements = {}
        result = self.formula.parseString(formula)
        self.parseresult(result, 1.)
        return self.elements.keys(), self.elements.values()


class CompoundFromFormula(compound.Compound):
    """Interface to a compound defined by a chemical formula
    """

    def __init__(self, formula, density=None, name=None):
        """
        Args:
            formula(str): chemical formula
            density(float): compound density (g/cm^3) 
            name(Optional[str]): compound name
        """
        p = FormulaParser()
        elements, mults = p.eval(formula)
        if not elements:
            raise ValueError(
                "Chemical formula {} is not valid".format(formula))
        if name is None:
            name = formula
        super(CompoundFromFormula, self).__init__(
            elements, mults, types.fraction.mole, density, name=name)


def factory(name):
    bsingle = instance.isstring(name)
    if bsingle:
        name = [name]
    try:
        ret = [CompoundFromFormula(s) for s in name]
        if bsingle:
            return ret[0]
        else:
            return ret
    except:
        return None


def search(name):
    return name
