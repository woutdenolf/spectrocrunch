# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.tri as plttri
import mpl_toolkits.axisartist as pltaa
import matplotlib.colors as pltcolors

from . import ternary_diagram
from . import colormap


def ChromaticityTriangle(fig, rect, ternaryinfo, debug=False, additive=True):

    ax = pltaa.Axes(fig, rect, frameon=False)
    if not debug:
        ax.axis["top", "right", "left", "bottom"].set_visible(False)
    H = np.sqrt(3) / 2
    ax.set_xlim(0, 1)
    ax.set_ylim(0, H)

    # Ternary grid
    n = 50
    pright = np.linspace(0, 1, n)
    ptop = np.linspace(0, 1, n)
    pright, ptop = np.meshgrid(pright, ptop)
    pright = pright.flatten()
    ptop = ptop.flatten()
    s = pright + ptop
    pleft = 1 - s

    r = 0.5 / n
    ind = (pleft >= 0) & ((s > (1 - r)) | (s < (1 + r)))
    pright = pright[ind]
    ptop = ptop[ind]
    pleft = pleft[ind]

    x = pright + ptop / 2.0
    y = H * ptop

    if additive:
        M = np.maximum.reduce([ptop, pleft, pright])
        ptop /= M
        pleft /= M
        pright /= M

    # Grid colors
    z = (
        np.round(pright * 255).astype(int) << ternaryinfo.right.shift
        | np.round(ptop * 255).astype(int) << ternaryinfo.top.shift
        | np.round(pleft * 255).astype(int) << ternaryinfo.left.shift
    )

    # Ternary grid
    triangles = plttri.Triangulation(x, y)
    ax.tripcolor(
        triangles,
        z,
        shading="gouraud",
        cmap=colormap.RGBcolormap(),
        norm=pltcolors.Normalize(vmin=0.0, vmax=2.0 ** 24 - 1),
    )

    return ax
