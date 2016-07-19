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

execfile("initcctbx.py")


# Don't use the installed version
import os, sys
sys.path.insert(1,os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from spectrocrunch.materials.compoundfromformula import compoundfromformula as compound
from spectrocrunch.materials.mixture import mixture
from spectrocrunch.materials.types import fractionType
import numpy as np

if __name__ == '__main__':

    compound1 = compound("La2O3",0,name="La2O3")
    compound2 = compound("SrO",0,name="SrO")
    compound3 = compound("Co2O3",0,name="Co2O3")
    compound4 = compound("Fe2O3",0,name="Fe2O3")

    m = mixture([compound1,compound2,compound3,compound4],[1,1,1,1],fractionType.weight)

    print(compound4.weightfractions())
    print("")
    elements = m.elemental_weightfractions()
    print("")
    for e in elements:
        print(e,e.MM,elements[e])

    

