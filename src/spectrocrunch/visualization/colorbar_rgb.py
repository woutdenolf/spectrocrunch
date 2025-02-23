import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as pltcolors

from . import ternary_diagram
from . import chromaticity_triangle


def triangle(
    fig=None,
    vmin=[0, 0, 0],
    vmax=[1, 1, 1],
    names=["R", "G", "B"],
    rect=[0.1, 0.1, 0.8, 0.8],
    grid=True,
    compositions=None,
):

    names = ["None" if name is None else name for name in names]

    right = ternary_diagram.TernaryComponent(
        name=names[0], min=vmin[0], max=vmax[0], color="r", shift=0
    )
    top = ternary_diagram.TernaryComponent(
        name=names[1], min=vmin[1], max=vmax[1], color="g", shift=8
    )
    left = ternary_diagram.TernaryComponent(
        name=names[2], min=vmin[2], max=vmax[2], color="b", shift=16
    )

    ternaryinfo = ternary_diagram.TernaryInfo(left=left, right=right, top=top)

    if fig is None:
        fig = plt.figure()

    items = {}
    items["colorbar"] = chromaticity_triangle.ChromaticityTriangle(
        fig, rect, ternaryinfo, additive=True
    )
    items["colorbar_axleft"] = ternary_diagram.axesLeft(fig, rect, ternaryinfo)
    items["colorbar_axtop"] = ternary_diagram.axesTop(fig, rect, ternaryinfo)
    items["colorbar_axright"] = ternary_diagram.axesRight(fig, rect, ternaryinfo)

    for item in items.values():
        fig.add_axes(item)

    if grid:
        ternary_diagram.TernaryGrid(items["colorbar"], ternaryinfo)

    if compositions:
        names, compositions, colors = zip(*compositions)
        right, top, left = zip(*compositions)
        ternary_diagram.TernaryLegend(
            items["colorbar"],
            ternaryinfo,
            np.array(left),
            np.array(right),
            np.array(top),
            names,
            colors,
        )

    return items


def bars(
    ax=None,
    vmin=[0, 0, 0],
    vmax=[1, 1, 1],
    names=["R", "G", "B"],
    norms=[None, None, None],
):
    cms = ([(0, 0, 0), (1, 0, 0)], [(0, 0, 0), (0, 1, 0)], [(0, 0, 0), (0, 0, 1)])

    pad = 0.05
    ret = {}

    for mi, ma, name, norm, cm in zip(vmin, vmax, names, norms, cms):
        if norm is None:
            norm = pltcolors.Normalize(vmin=mi, vmax=ma)
        cm = pltcolors.LinearSegmentedColormap.from_list(name, cm)
        sm = plt.cm.ScalarMappable(cmap=cm, norm=norm)
        sm._A = []
        cb = plt.colorbar(sm, ax=ax, fraction=0.15, pad=pad)
        label = cb.ax.text(0, 1.05, name, transform=cb.ax.transAxes)
        ret["colorbar{}".format(name)] = cb
        ret["colorbar{}_title".format(name)] = label
        pad += 0.05

    return ret
