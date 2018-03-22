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
import matplotlib.tri as plttri
import mpl_toolkits.axisartist as pltaa
import matplotlib.colors as pltcolors

from . import ternary_diagram
from . import colormap

def ChromaticityTriangle(fig,rect,ternaryinfo,debug=False,additive=True):

    ax = pltaa.Axes(fig,rect,frameon=False)
    if not debug:
        ax.axis["top", "right","left", "bottom"].set_visible(False)
    H = np.sqrt(3)/2
    ax.set_xlim(0,1)
    ax.set_ylim(0,H)
    
    # Ternary grid
    n = 50
    pright = np.linspace(0,1,n)
    ptop = np.linspace(0,1,n)
    pright,ptop = np.meshgrid(pright,ptop)
    pright = pright.flatten()
    ptop = ptop.flatten()
    s = pright+ptop
    pleft = 1-s
    
    r = 0.5/n
    ind = (pleft>=0) & ((s>(1-r)) | (s<(1+r)))
    pright = pright[ind]
    ptop = ptop[ind]
    pleft = pleft[ind]
    
    x = pright + ptop/2.
    y = H*ptop
    
    if additive:
        M = np.maximum.reduce([ptop,pleft,pright])
        ptop /= M
        pleft /= M
        pright /= M

    # Grid colors
    z = np.round(pright*255).astype(int) << ternaryinfo.right.shift |\
         np.round(ptop*255).astype(int) << ternaryinfo.top.shift  |\
         np.round(pleft*255).astype(int) << ternaryinfo.left.shift 

    # Ternary grid
    triangles = plttri.Triangulation(x, y)
    ax.tripcolor(triangles, z, shading='gouraud', cmap=colormap.RGBcolormap(),\
                        norm=pltcolors.Normalize(vmin=0.,vmax=2.**24-1))

    if debug:
        ax.triplot(triangles, 'bo-', lw=1, color=(1,0,0))
        ternary_diagram.TernaryPoint(ax,ternaryinfo,ternary_diagram.TernaryCoordinates(left=0.2,right=0.6,top=0.2))

    return ax
    
