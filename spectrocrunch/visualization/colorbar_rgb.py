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

import numpy as np
import matplotlib.pyplot as plt

from . import ternary_diagram
from . import chromaticity_triangle

def colorbar_rgb(fig=None,vmin=[0,0,0],vmax=[1,1,1],names=['red','green','blue'],rect=[0.1, 0.1, 0.8, 0.8],grid=True):

    names = ["None" if name is None else name for name in names]

    right = ternary_diagram.TernaryComponent(name=names[0],min=vmin[0],max=vmax[0],color='r',shift=0)
    top = ternary_diagram.TernaryComponent(name=names[1],min=vmin[1],max=vmax[1],color='g',shift=8)
    left = ternary_diagram.TernaryComponent(name=names[2],min=vmin[2],max=vmax[2],color='b',shift=16)
    
    ternaryinfo = ternary_diagram.TernaryInfo(left=left,right=right,top=top)

    if fig is None:
        fig = plt.figure()
        
    items = []
    items.append(chromaticity_triangle.ChromaticityTriangle(fig,rect,ternaryinfo,additive=True))
    items.append(ternary_diagram.axesLeft(fig,rect,ternaryinfo))
    items.append(ternary_diagram.axesTop(fig,rect,ternaryinfo))
    items.append(ternary_diagram.axesRight(fig,rect,ternaryinfo))

    for item in items:
        fig.add_axes(item)
    
    if grid:
        ternary_diagram.TernaryGrid(items[0],ternaryinfo)

    return fig,items
    
