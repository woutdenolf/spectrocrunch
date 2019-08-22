# -*- coding: utf-8 -*-

import numpy as np
import collections
import matplotlib.colors as pltcolors
import matplotlib.transforms as plttransforms
import matplotlib.patches as pltpatches
import mpl_toolkits.axisartist as plta
import mpl_toolkits.axisartist.grid_helper_curvelinear as pltahelper

CartesianCoordinates = collections.namedtuple(
    'CartesianCoordinates', ['x', 'y'])
TernaryCoordinates = collections.namedtuple(
    'TernaryCoordinates', ['left', 'right', 'top'])
TernaryComponent = collections.namedtuple(
    'TernaryComponent', ['name', 'min', 'max', 'color', 'shift'])
TernaryInfo = collections.namedtuple('TernaryInfo', ['left', 'right', 'top'])


def transformations(M):
    Mi = M.inverted()

    def transform(m, x, y):
        x = np.atleast_1d(x)
        y = np.atleast_1d(y)
        xy = np.stack((x, y), axis=1)
        xy = m.transform(xy)
        return xy[:, 0], xy[:, 1]

    def f(x, y): return transform(M, x, y)
    def fi(x, y): return transform(Mi, x, y)
    return f, fi


def axesLeft(fig, rect, ternaryinfo, grid=False, debug=False):
    xname = ternaryinfo.left.name
    vminx, vmaxx = ternaryinfo.left.min, ternaryinfo.left.max
    yname = ternaryinfo.top.name
    vminy, vmaxy = ternaryinfo.top.max, ternaryinfo.top.min

    ox = vminx
    oy = vminy
    sx = float(vmaxx-vminx)
    sy = float(vmaxy-vminy)

    M = plttransforms.Affine2D()
    M += plttransforms.Affine2D().translate(-ox, -oy)
    M += plttransforms.Affine2D().scale(1, sx/sy)
    M += plttransforms.Affine2D().skew_deg(15, 15)
    M += plttransforms.Affine2D().rotate_deg(-120-15)
    f, fi = transformations(M)

    grid_helper = pltahelper.GridHelperCurveLinear((f, fi))
    ax = plta.Axes(fig, rect, grid_helper=grid_helper, frameon=False)
    ax.axis["top", "right", "left", "bottom"].set_visible(False)

    ax.axis["X"] = ax.new_floating_axis(1, oy)
    ax.axis["X"].label.set_text(xname)
    ax.axis["X"].label.set_color(ternaryinfo.left.color)

    xmin1, ymin1 = f(vmaxx, vminy)
    xmin1, ymin1 = xmin1[0], ymin1[0]
    xmax1, _ = f(vminx, vmaxy)
    xmax1, ymax1 = xmax1[0], 0
    ax.set_xlim(xmin1, xmax1)
    ax.set_ylim(ymin1, ymax1)

    if debug:
        ax.axis["Y"] = ax.new_floating_axis(0, ox)
        ax.axis["Y"].label.set_text(yname)
        ax.axis["Y"].toggle(ticks=False)

        ax.plot(*f([vminx, vmaxx], [vminy, vminy]), marker="o", markersize=10)
        ax.plot(*f([vminx, vminx], [vminy, vmaxy]), marker="+", markersize=10)

    if grid:
        ax.grid(True, axis="x", zorder=0, color=ternaryinfo.left.color)

    return ax


def axesTop(fig, rect, ternaryinfo, grid=False, debug=False):
    xname = ternaryinfo.right.name
    vminx, vmaxx = ternaryinfo.right.max, ternaryinfo.right.min
    yname = ternaryinfo.top.name
    vminy, vmaxy = ternaryinfo.top.min, ternaryinfo.top.max

    ox = vminx
    oy = vminy
    sx = float(vmaxx-vminx)
    sy = float(vmaxy-vminy)

    M = plttransforms.Affine2D()
    M += plttransforms.Affine2D().translate(-ox, -oy)
    M += plttransforms.Affine2D().scale(1, sx/sy)
    M += plttransforms.Affine2D().skew_deg(15, 15)
    M += plttransforms.Affine2D().rotate_deg(-15)
    f, fi = transformations(M)

    grid_helper = pltahelper.GridHelperCurveLinear((f, fi))
    ax = plta.Axes(fig, rect, grid_helper=grid_helper, frameon=False)
    ax.axis["top", "right", "left", "bottom"].set_visible(False)

    ax.axis["Y"] = ax.new_floating_axis(0, ox)
    ax.axis["Y"].label.set_text(yname)
    ax.axis["Y"].label.set_color(ternaryinfo.top.color)

    xmax1, ymin1 = f(vminx, vminy)
    xmax1, ymin1 = xmax1[0], ymin1[0]
    xmin1, _ = f(vmaxx, vminy)
    xmin1 = xmin1[0]
    _, ymax1 = f(vminx, vmaxy)
    ymax1 = ymax1[0]
    ax.set_xlim(xmin1, xmax1)
    ax.set_ylim(ymin1, ymax1)

    if debug:
        ax.axis["X"] = ax.new_floating_axis(1, oy)
        ax.axis["X"].label.set_text(xname)

        ax.plot(*f([vminx, vmaxx], [vminy, vminy]), marker="o", markersize=10)
        ax.plot(*f([vminx, vminx], [vminy, vmaxy]), marker="+", markersize=10)

    if grid:
        ax.grid(True, axis="y", zorder=0, color=ternaryinfo.top.color)

    return ax


def axesRight(fig, rect, ternaryinfo, grid=False, debug=False):
    xname = ternaryinfo.right.name
    vminx, vmaxx = ternaryinfo.right.max, ternaryinfo.right.min
    yname = ternaryinfo.top.name
    vminy, vmaxy = ternaryinfo.top.min, ternaryinfo.top.max

    ox = vminx
    oy = vminy
    sx = float(vmaxx-vminx)
    sy = float(vmaxy-vminy)

    M = plttransforms.Affine2D()
    M += plttransforms.Affine2D().translate(-ox, -oy)
    M += plttransforms.Affine2D().scale(1, sx/sy)
    M += plttransforms.Affine2D().skew_deg(-15, -15)
    M += plttransforms.Affine2D().rotate_deg(15)
    f, fi = transformations(M)

    grid_helper = pltahelper.GridHelperCurveLinear((f, fi))
    ax = plta.Axes(fig, rect, grid_helper=grid_helper, frameon=False)
    ax.axis["top", "right", "left", "bottom"].set_visible(False)

    ax.axis["X"] = ax.new_floating_axis(1, oy)
    ax.axis["X"].label.set_text(xname)
    ax.axis["X"].label.set_color(ternaryinfo.right.color)

    xmax1, ymin1 = f(vminx, vminy)
    xmax1, ymin1 = xmax1[0], ymin1[0]
    xmin1, _ = f(vmaxx, vminy)
    xmin1 = xmin1[0]
    _, ymax1 = f(vminx, vmaxy)
    ymax1 = ymax1[0]
    ax.set_xlim(xmin1, xmax1)
    ax.set_ylim(ymin1, ymax1)

    if debug:
        ax.axis["Y"] = ax.new_floating_axis(0, ox)
        ax.axis["Y"].label.set_text(yname)

        ax.plot(*f([vminx, vmaxx], [vminy, vminy]), marker="o", markersize=10)
        ax.plot(*f([vminx, vminx], [vminy, vmaxy]), marker="+", markersize=10)

    if grid:
        ax.grid(True, axis="x", zorder=0, color=ternaryinfo.right.color)

    return ax


def NormTernary(p, ternaryinfo):
    def clip(x): return np.clip(x, 0, 1)

    mi, ma = ternaryinfo.top.min, ternaryinfo.top.max
    top = clip((p.top-mi)/(ma-mi))

    mi, ma = ternaryinfo.left.min, ternaryinfo.left.max
    left = clip((p.left-mi)/(ma-mi))

    mi, ma = ternaryinfo.right.min, ternaryinfo.right.max
    right = clip((p.right-mi)/(ma-mi))

    return TernaryCoordinates(left=left, right=right, top=top)


def TernaryToCartesian(p):
    H = np.sqrt(3)/2
    tot = p.left+p.right+p.top
    return CartesianCoordinates(x=(2*p.right+p.top)/(2.*tot), y=H*p.top/tot)


def CartesianToTernary(p):
    H = np.sqrt(3)/2
    top = p.y/H
    right = p.x-top/2
    left = 1-top-right
    return TernaryCoordinates(left=left, right=right, top=top)


def TernaryPoint(ax, ternaryinfo, point):
    H = np.sqrt(3)/2

    x0, y0 = TernaryToCartesian(point)

    x1, y1 = TernaryToCartesian(
        point._replace(left=point.left+point.top, top=0))
    ax.plot([x0, x1], [y0, y1], color=ternaryinfo.right.color)

    x1, y1 = TernaryToCartesian(point._replace(
        top=point.top+point.right, right=0))
    ax.plot([x0, x1], [y0, y1], color=ternaryinfo.left.color)

    x1, y1 = TernaryToCartesian(point._replace(
        right=point.right+point.left, left=0))
    ax.plot([x0, x1], [y0, y1], color=ternaryinfo.top.color)

    ax.plot([x0], [y0], marker='o', markersize=10, color="#000000")


def TernaryGrid(ax, ternaryinfo, n=5):

    p = np.repeat(np.linspace(0, 1, n+2)[1:-1], 2).reshape((n, 2))

    right = p.copy()
    left = 1-right
    top = 1-right
    left[:, 0] = 0
    top[:, 1] = 0
    pts = TernaryCoordinates(left=left, right=right, top=top)
    x, y = TernaryToCartesian(pts)
    for xi, yi in zip(x, y):
        ax.plot(xi, yi, color=ternaryinfo.right.color)

    left = p.copy()
    right = 1-left
    top = 1-left
    right[:, 0] = 0
    top[:, 1] = 0
    pts = TernaryCoordinates(left=left, right=right, top=top)
    x, y = TernaryToCartesian(pts)
    for xi, yi in zip(x, y):
        ax.plot(xi, yi, color=ternaryinfo.left.color)

    top = p.copy()
    right = 1-top
    left = 1-top
    right[:, 0] = 0
    left[:, 1] = 0
    pts = TernaryCoordinates(left=left, right=right, top=top)
    x, y = TernaryToCartesian(pts)
    for xi, yi in zip(x, y):
        ax.plot(xi, yi, color=ternaryinfo.top.color)


def TernaryLegend(ax, ternaryinfo, left, right, top, names, colors):
    pts = TernaryCoordinates(left=left, right=right, top=top)
    pts = NormTernary(pts, ternaryinfo)

    if False:
        x, y = TernaryToCartesian(pts)
        for xi, yi, name, c in zip(x, y, names, colors):
            ax.plot(x, y, 'o', color="#000000")
    else:
        TernaryPoint(ax, ternaryinfo, pts)

    patches = [pltpatches.Patch(color=color, label=label)
               for label, color in zip(names, colors)]
    ax.legend(patches, names, loc="upper right", frameon=False)
