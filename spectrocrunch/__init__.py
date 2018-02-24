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

try:
    from ._version import version as __version__
except ImportError:
    import os
    __version__ = "Local version ({})".format(os.path.dirname(os.path.abspath(__file__)))

##### Make sure there is always a logger #####
import logging
logging.basicConfig()

##### Initialize pint #####
import pint
ureg = pint.UnitRegistry()

import scipy.constants
r = scipy.constants.physical_constants["classical electron radius"]
ureg.re = ureg.Quantity(r[0],r[1])

##### Initialize xraylib #####
import xraylib
xraylib.XRayInit()
xraylib.SetErrorMessages(0)

# Code <-> Name: one-to-one
xraylib.code_to_shell = {xraylib.__dict__[s]:s.split('_')[0] for s in xraylib.__dict__.keys() if s.endswith("_SHELL")}
xraylib.shell_to_code = {v: k for k, v in xraylib.code_to_shell.items()}
xraylib.shell_min = min(xraylib.code_to_shell)
xraylib.shell_max = max(xraylib.code_to_shell)

# Code <-> Name: one-to-many
xraylib.line_to_code = {s.split('_')[0]:xraylib.__dict__[s] for s in xraylib.__dict__.keys() if s.endswith("_LINE")}
xraylib.code_to_line = {code:[name for name,code2 in xraylib.line_to_code.items() if code==code2] for code in set(xraylib.line_to_code.values())}
xraylib.line_min = min(xraylib.code_to_line)
xraylib.line_max = max(xraylib.code_to_line)
   
# Composites
xraylib.composites ={xraylib.KA_LINE:[xraylib.KL3_LINE, xraylib.KL2_LINE, xraylib.KL1_LINE],\
                      xraylib.KB_LINE:[xraylib.KP5_LINE, xraylib.KP4_LINE, xraylib.KP3_LINE, xraylib.KP2_LINE, xraylib.KP1_LINE]+\
                                      [xraylib.KO7_LINE, xraylib.KO6_LINE, xraylib.KO5_LINE, xraylib.KO4_LINE, xraylib.KO3_LINE, xraylib.KO2_LINE, xraylib.KO1_LINE]+\
                                      [xraylib.KN7_LINE, xraylib.KN6_LINE, xraylib.KN5_LINE, xraylib.KN4_LINE, xraylib.KN3_LINE, xraylib.KN2_LINE, xraylib.KN1_LINE]+\
                                      [xraylib.KM5_LINE, xraylib.KM4_LINE, xraylib.KM3_LINE, xraylib.KM2_LINE, xraylib.KM1_LINE],\
                      xraylib.KP_LINE:[xraylib.KP5_LINE, xraylib.KP4_LINE, xraylib.KP3_LINE, xraylib.KP2_LINE, xraylib.KP1_LINE],\
                      xraylib.KO_LINE:[xraylib.KO7_LINE, xraylib.KO6_LINE, xraylib.KO5_LINE, xraylib.KO4_LINE, xraylib.KO3_LINE, xraylib.KO2_LINE, xraylib.KO1_LINE],\
                      xraylib.LA_LINE:[xraylib.L3M5_LINE, xraylib.L3M4_LINE],\
                      xraylib.LB_LINE:[xraylib.L3O4_LINE, xraylib.L3O3_LINE, xraylib.L3N5_LINE, xraylib.L3N1_LINE, xraylib.L2M4_LINE, xraylib.L1M3_LINE, xraylib.L1M2_LINE],\
                      xraylib.L1N67_LINE:[xraylib.L1N7_LINE,xraylib.L1N6_LINE],\
                      xraylib.L1O45_LINE:[xraylib.L1O5_LINE,xraylib.L1O4_LINE],\
                      xraylib.L1P23_LINE:[xraylib.L1P3_LINE,xraylib.L1P2_LINE],\
                      xraylib.L2P23_LINE:[xraylib.L2P3_LINE,xraylib.L2P2_LINE],\
                      xraylib.L3O45_LINE:[xraylib.L3O5_LINE,xraylib.L3O4_LINE],\
                      xraylib.L3P23_LINE:[xraylib.L3P3_LINE,xraylib.L3P2_LINE],\
                      xraylib.L3P45_LINE:[xraylib.L3P5_LINE,xraylib.L3P4_LINE]\
                     }

xraylib.rcomposites = {}
for k,v in xraylib.composites.items():
    for l in v:
        if l in xraylib.rcomposites:
            xraylib.rcomposites[l].append(k)
            xraylib.rcomposites[l].sort()
        else:
            xraylib.rcomposites[l]=[k]

