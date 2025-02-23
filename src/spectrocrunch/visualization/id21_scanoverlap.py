from ..io.spec import spec
import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from ..math.utils import logscale
import warnings


def show(
    x,
    y,
    images,
    xp,
    yp,
    xlabel,
    ylabel,
    names,
    transpose=False,
    flipvert=False,
    fliphor=False,
    color="#ffffff",
    defaultorigin=False,
    printpos=False,
    outname=None,
):
    """
    Args:
        x(np.array): horizontal coordinates
        y(np.array): vertical coordinates
        images(np.array): image
        xp(np.array): marker horizontal coord.
        yp(np.array): marker vertical coord.
        xlabel(str):
        ylabel(str):
        names(list(str)):
    """

    # Make monotonically increasing (required by interp2d)
    ind = np.argsort(x)
    x = x[ind]
    images = images[:, :, ind]
    ind = np.argsort(y)
    y = y[ind]
    images = images[:, ind, :]
    nimg = images.shape[0]

    # New grid
    xnew = np.linspace(x[0], x[-1], len(x))
    ynew = np.linspace(y[0], y[-1], len(y))

    # Interpolate
    for i in range(nimg):
        # Use another algorithm
        f = interpolate.interp2d(x, y, images[i, ...], kind="cubic")
        images[i, ...] = np.clip(f(xnew, ynew), 0, 1)

    # Plot range
    dx = (xnew[1] - xnew[0]) / 2.0
    dy = (ynew[1] - ynew[0]) / 2.0
    extent = (x[0] - dx, x[-1] + dx, y[0] - dy, y[-1] + dy)
    origin = "lower"

    # Transpose
    if transpose:
        extent = (extent[2], extent[3], extent[0], extent[1])
        images = images.transpose((0, 2, 1))
        xp, yp = yp, xp
        xlabel, ylabel = ylabel, xlabel

    # Flip vertical
    if flipvert:
        extent = (extent[0], extent[1], extent[3], extent[2])
        images = images[:, ::-1, :]

    # Flip horizontal
    if fliphor:
        extent = (extent[1], extent[0], extent[2], extent[3])
        images = images[:, :, ::-1]

    # Origin left bottom
    if defaultorigin:
        ind = [0, 1, 2, 3]
        if extent[1] < extent[0]:
            # extent[0] ... xp  .......... extent[1]
            # extent[1] ... xp  .......... extent[0]
            xp = extent[1] + extent[0] - xp
            ind[0] = 1
            ind[1] = 0
        if extent[3] < extent[2]:
            ind[2] = 3
            ind[3] = 2
            yp = extent[3] + extent[2] - yp
        extent = (extent[ind[0]], extent[ind[1]], extent[ind[2]], extent[ind[3]])

    # Show
    if printpos:
        print(extent)
        print(np.vstack((xp, yp)).T)

    # RGB for plotting
    if transpose:
        rgb = np.zeros((len(xnew), len(ynew), 3))
    else:
        rgb = np.zeros((len(ynew), len(xnew), 3))
    for i in range(3):
        rgb[..., i] = images[i, ...]
    # rgb = images[0:3,...].transpose((1,2,0))

    # Plot
    plt.figure(1)
    plt.clf()
    _ = plt.imshow(
        rgb, extent=extent, origin=origin, interpolation="nearest", aspect=1
    )  # ,cmap=plt.get_cmap("gray")
    axes = plt.gca()
    axes.set_xlabel(xlabel)
    axes.set_ylabel(ylabel)
    xlim, ylim = axes.get_xlim(), axes.get_ylim()

    fontsize = 12
    s = fontsize / 2
    axes.scatter(xp, yp, marker="o", s=s, color=color)
    for i in range(len(names)):
        # try:
        #    rgbi = rgb[int(np.round(xp[i])),int(np.round(yp[i])),:]*255

        # print(rgbi[0]*0.299 + rgbi[1]*0.587 + rgbi[2]*0.114)
        #    if (rgbi[0]*0.299 + rgbi[1]*0.587 + rgbi[2]*0.114) > 100:
        #        color = '#000000'
        #    else:
        #        color = '#ffffff'
        # except Exception:
        #    color = '#ffffff'

        # color = '#%02x%02x%02x' % tuple(255-rgbi)
        # axes.scatter(xp[i], yp[i], marker='o',s=s,color = color)

        if names[i] is not None:
            axes.annotate(
                names[i], xy=(xp[i], yp[i]), xytext=(xp[i] + dx, yp[i]), color=color
            )

    axes.set_xlim(xlim)
    axes.set_ylim(ylim)
    if outname is None:
        plt.show()
    else:
        plt.savefig(outname, bbox_inches="tight", dpi=300)


def plot(
    hdf5filename,
    grps,
    specfilename,
    specnumbers,
    offsamy,
    offsamz,
    transpose=False,
    flipvert=True,
    fliphor=False,
    defaultorigin=False,
    showlabels=False,
    color="#ffffff",
    printpos=False,
    outname=None,
    log=False,
):
    """
    Args:
        hdf5filename(str)
        grps(dict): keys must be 0, 1 or 2 (r, g, b)
        specfilename(str)
        specnumbers(list(int))
        offhor(float)
        offvert(float)
    """

    oh5 = h5py.File(hdf5filename)

    # Prepare global coordinates
    dim1off = 0.0
    dim1name = "samz"
    dim1mult = 1
    dim2off = 0.0
    dim2name = "samy"
    dim2mult = 1
    try:
        ocoord = oh5["stackinfo"]
    except KeyError:
        warnings.warn(
            '"coordinates" is deprecated and should be replaced by "stackinfo"',
            DeprecationWarning,
        )
        ocoord = oh5["coordinates"]
    for f in ocoord:
        if f == "samz":
            dim1off = ocoord[f].value * 1000
            dim1name = "sampz"
            dim1mult = 1
        if f == "sampz":
            dim1off = ocoord[f].value
            dim1name = "samz"
            dim1mult = 1000
        if f == "samy":
            dim2off = ocoord[f].value * 1000
            dim2name = "sampy"
            dim2mult = 1
        if f == "sampy":
            dim2off = ocoord[f].value
            dim2name = "samy"
            dim2mult = 1000

    # Get image with axes in micron
    for i in grps:
        ogrp = oh5[grps[i]["path"]]
        odset = ogrp[ogrp.attrs["signal"]]
        dim1 = dim1off[grps[i]["ind"]] + ogrp[dim1name].value * dim1mult
        dim2 = dim2off[grps[i]["ind"]] + ogrp[dim2name].value * dim2mult
        idim1 = ogrp.attrs[dim1name + "_indices"]
        idim2 = ogrp.attrs[dim2name + "_indices"]
        if idim2 != 0 and idim1 != 0:
            img = odset[grps[i]["ind"], ...]
        elif idim2 != 1 and idim1 != 1:
            img = odset[:, grps[i]["ind"], :]
        else:
            img = odset[..., grps[i]["ind"]]
        img[np.isnan(img)] = np.nanmin(img)
        if idim1 > idim2:
            img = img.T
        if i == 0:
            images = np.zeros((3,) + img.shape, dtype=img.dtype)

        if log:
            img = logscale(img)

        mi = np.min(img)
        ma = np.max(img)
        d = ma - mi
        mi += d * grps[i]["lo"]
        ma -= d * (1 - grps[i]["hi"])
        img -= mi
        img /= ma
        img = np.clip(img, 0, 1)

        images[i, ...] = img

    oh5.close()

    # XANES positions
    ospec = spec(specfilename)
    motors = ["samz", "sampz", "samy", "sampy"]
    n = len(specnumbers)
    pdim1 = np.empty(n)
    pdim2 = np.empty(n)
    if not hasattr(offsamz, "__len__"):
        offsamz = [offsamz] * n
    if not hasattr(offsamy, "__len__"):
        offsamy = [offsamy] * n
    for i in range(n):
        v = ospec.getmotorvalues(specnumbers[i], motors)
        if printpos:
            print("Spec number {}".format(i))
            for a, b in zip(motors, v):
                print(" {} = {}".format(a, b))
        pdim1[i] = v[0] * 1000 + v[1] + offsamz[i]
        pdim2[i] = v[2] * 1000 + v[3] + offsamy[i]

    # Make axes values readable
    m1 = min(dim1)
    m2 = min(dim2)
    dim1 -= m1
    dim2 -= m2
    pdim1 -= m1
    pdim2 -= m2

    # Plot
    if showlabels:
        names = [str(i) for i in specnumbers]
    else:
        names = [None] * n
    if defaultorigin:
        dim2label = r"x ($\mu$m)"
        dim1label = r"y ($\mu$m)"
    else:
        dim2label = r"y ($\mu$m)"
        dim1label = r"z ($\mu$m)"

    show(
        dim2,
        dim1,
        images,
        pdim2,
        pdim1,
        dim2label,
        dim1label,
        names,
        transpose=transpose,
        flipvert=flipvert,
        fliphor=fliphor,
        color=color,
        defaultorigin=defaultorigin,
        printpos=printpos,
        outname=outname,
    )
