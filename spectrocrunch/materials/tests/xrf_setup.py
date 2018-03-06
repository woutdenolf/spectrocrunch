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

from ...geometries import xrf as xrfgeometries
from ...sources import xray as xraysources
from ...detectors import xrf as xrfdetectors
from .. import compoundfromname
from .. import compoundfromformula
from .. import element
from .. import mixture
from .. import types
from .. import multilayer
from .. import pymca
from .. import xrfstandards

source = xraysources.factory("synchrotron")
detector = xrfdetectors.factory("leia")
geometry = xrfgeometries.factory("sxm120",detector=detector,
                                 source=source,detectorposition=-15.)

# Cover elements, compounds and mixtures
# Not too many lines for speed
hematite = compoundfromname.compoundfromname("hematite")
goethite = compoundfromname.compoundfromname("goethite")
mix = mixture.Mixture([goethite,hematite],[0.5,0.5],types.fractionType.weight,name="iron oxides")
calcite = compoundfromname.compoundfromname("calcite")
ca = element.Element("Ca")
sample = multilayer.Multilayer([ca,mix,calcite],[2e-5,7e-5,10e-5],geometry=geometry)

# Cover L and M lines
compound1 = compoundfromformula.CompoundFromFormula("PbCe",density=6.)
sample_complex = multilayer.Multilayer([ca,mix,compound1,calcite],[2e-5,3e-5,1e-5,10e-5],geometry=geometry)



