# -*- coding: utf-8 -*-

import numpy as np

# import scipy.ndimage
from scipy.stats import rv_continuous

from ..types import transformationType


def makeGaussian(x, y, x0, y0, sx, sy, rho, A):
    num = (
        (x - x0) ** 2 / sx**2
        - 2 * rho / (sx * sy) * (x - x0) * (y - y0)
        + (y - y0) ** 2 / sy**2
    )
    denom = 2 * (1 - rho**2)
    return A * np.exp(-num / denom)  # /(2*np.pi*sx*sy*np.sqrt(1-rho**2))


def makeGaussianAngle(x, y, x0, y0, sx, sy, angle, A):
    cosa = np.cos(angle)
    sina = np.sin(angle)
    sin2a = np.sin(2 * angle)
    varx = sx**2
    vary = sy**2
    a = cosa**2 / (2 * varx) + sina**2 / (2 * vary)
    b = -sin2a / (4 * varx) + sin2a / (4 * vary)
    c = sina**2 / (2 * varx) + cosa**2 / (2 * vary)
    num = a * (x - x0) ** 2 - 2 * b * (x - x0) * (y - y0) + c * (y - y0) ** 2
    return A * np.exp(-num)


def gettransformedimage(x, y, data, angle=False):
    ret = np.zeros(x.shape, dtype=float)
    if angle:
        proc = makeGaussianAngle
    else:
        proc = makeGaussian
    for x0, y0, sx, sy, rho, A in data:
        ret += proc(x, y, x0, y0, sx, sy, rho, A)
    # ret /= np.max(ret)
    return ret


# def getlargeimage(npeaks,nsigma,shape,subshape):
#    n = np.round(npeaks*min(shape)/min(subshape)).astype(int)
#    image = np.zeros(shape,dtype=np.float32)
#    np.random.seed(1)
#    x = shape[0]*np.random.random(n**2)
#    y = shape[1]*np.random.random(n**2)
#    image[x.astype(int), y.astype(int)] = 1
#    image = scipy.ndimage.gaussian_filter(image, sigma=min(subshape)/(npeaks*nsigma))
#    image /= np.max(image)


def makeGaussian1D(x, x0, sx, A):
    return A * np.exp(-((x - x0) ** 2) / (2.0 * sx**2))  # /(np.sqrt(2*np.pi)*sx)


def gettransformedvector(x, data, angle=False):
    ret = np.zeros(x.shape, dtype=float)
    for x0, sx, A in data:
        ret += makeGaussian1D(x, x0, sx, A)
    # ret /= np.max(ret)
    return ret


def translation(dx, dy):
    Mcoord = np.identity(3)
    Mcof = np.identity(3)
    Mcoord[0:2, 2] = [dx, dy]
    Mcof[0:2, 2] = [-dx, -dy]
    return Mcoord, Mcof


def rigid(a, tx, ty):
    if a < 0 or a > 1:
        raise ValueError("A rigid transformation is expected.")
    Mcoord = np.identity(3)
    Mcof = np.identity(3)
    b = np.sqrt(1 - a * a)
    Mcoord[0, 0:3] = [a, -b, tx]
    Mcoord[1, 0:3] = [b, a, ty]
    Mcof[0, 0:3] = [a, b, -a * tx - b * ty]
    Mcof[1, 0:3] = [-b, a, b * tx - a * ty]
    return Mcoord, Mcof


def similarity(a, b, tx, ty):
    Mcoord = np.identity(3)
    Mcof = np.identity(3)
    Mcoord[0, 0:3] = [a, -b, tx]
    Mcoord[1, 0:3] = [b, a, ty]
    s = float(a * a + b * b)
    Mcof[0, 0:3] = [a / s, b / s, -(a * tx + b * ty) / s]
    Mcof[1, 0:3] = [-b / s, a / s, (b * tx - a * ty) / s]
    return Mcoord, Mcof


def affine(a, b, c, d, tx, ty):
    Mcoord = np.identity(3)
    Mcof = np.identity(3)
    Mcoord[0, 0:3] = [a, b, tx]
    Mcoord[1, 0:3] = [c, d, ty]
    det = float(a * d - b * c)
    Mcof[0, 0:3] = [d / det, -b / det, (b * ty - d * tx) / det]
    Mcof[1, 0:3] = [-c / det, a / det, (c * tx - a * ty) / det]
    return Mcoord, Mcof


def projective(a, b, c, d, tx, ty, px, py):
    Mcoord = np.identity(3)
    Mcof = np.identity(3)
    Mcoord[0, 0:3] = [a, b, tx]
    Mcoord[1, 0:3] = [c, d, ty]
    Mcoord[2, 0:3] = [px, py, 1]
    det = float(a * d - a * py * ty - b * c + b * px * ty + c * py * tx - d * px * tx)
    Mcof[0, 0:3] = [d - py * ty, py * tx - b, b * ty - d * tx]
    Mcof[1, 0:3] = [px * ty - c, a - px * tx, c * tx - a * ty]
    Mcof[2, 0:3] = [c * py - d * px, b * px - a * py, a * d - b * c]
    Mcof /= det
    return Mcoord, Mcof


def transformation(t, n, subpixel=True):
    if t == transformationType.rigid:
        # total rotation is 40 degrees
        a = np.around(np.cos(40.0 / 180.0 * np.pi / (n - 1)), 3)
        # a = np.sqrt(3)/2 # between -1 and 1
    elif t == transformationType.similarity:
        a = np.around(np.cos(40.0 / 180.0 * np.pi / (n - 1)), 3)
        b = np.around(np.sin(40.0 / 180.0 * np.pi / (n - 1)), 3)
        a *= 1.1
        b *= 1.1
    else:
        a = np.around(np.cos(40.0 / 180.0 * np.pi / (n - 1)), 3)
        b = np.around(np.sin(40.0 / 180.0 * np.pi / (n - 1)), 3)
        a *= 1.1
        b *= 1.2
        c = -b * 0.9
        d = a * 0.9

    if subpixel:
        tx = 1.5
        ty = -1.7
    else:
        tx = 1.0
        ty = -4.0  # TODO; -2
    px = 0.001
    py = -0.001

    if t == transformationType.translation:
        return translation(tx, ty)
    elif t == transformationType.rigid:
        return rigid(a, tx, ty)
    elif t == transformationType.similarity:
        return similarity(a, b, tx, ty)
    elif t == transformationType.affine:
        return affine(a, b, c, d, tx, ty)
    elif t == transformationType.projective:
        return projective(a, b, c, d, tx, ty, px, py)
    else:
        raise NotImplementedError("This transformation to has not been implemented.")


class rvgen(rv_continuous):
    def _pdf(self, x):
        H = 1 / np.sqrt(2.0 * np.pi)
        return H * np.exp(-(x**2) / 2.0)

        return H * (1 - np.exp(-(x**2) / 2.0))


def random(a, b, n):
    return a + (b - a) * np.random.random(n)


def genpos1d(a, b, npeaks):
    x = random(a, b, npeaks)
    dr = (b - a) * 0.1
    xm = (a + b) / 2.0
    ind = abs(x - xm) > dr
    x = x[ind]
    npeaks = sum(ind)
    return x, npeaks


def genpos2d(a, b, c, d, npeaks):
    x = random(a, b, npeaks)
    y = random(c, d, npeaks)
    dr = max(b - a, d - c) * 0.1
    xm = (a + b) / 2.0
    ym = (c + d) / 2.0
    r = ((x - xm) ** 2 + (y - ym) ** 2) ** 0.5
    ind = r > dr
    x = x[ind]
    y = y[ind]
    npeaks = sum(ind)
    return x, y, npeaks


def data(
    transfotype,
    ndim1=61,  # TODO 61
    ndim2=57,
    nimages=4,
    nstacks=2,
    ndim=100,
    vector=False,
    transposed=False,
    realistic=True,
    subpixel=True,
    plot=False,
    inverse=False,
):
    """
    Returns:
        3-tuple: list of image stacks, change-of-coordinate matrix between the subsequent images, stack dimensions
    """
    if vector and transfotype != transformationType.translation:
        raise ValueError("Vectors can only be shifted.")

    # Transformation between images
    Mcoord, Mcof = transformation(transfotype, nimages, subpixel=subpixel)
    stackdim = 0

    if vector:
        # Shape of a subimage
        if transposed:
            dimi = 0
            subshape = (ndim, 1)
        else:
            dimi = 1
            subshape = (1, ndim)
        Mcoord[dimi, 2] = 0
        Mcof[dimi, 2] = 0
        subxmin = -subshape[dimi] // 2
        subxmax = subshape[dimi] // 2
        if transposed:
            subshape = (subxmax - subxmin + 1, 1)
        else:
            subshape = (1, subxmax - subxmin + 1)

        # Determine shape of large image (transform corners)
        d = Mcof[1 - dimi, 2] * (nimages - 1)
        xmin = int(min(subxmin, subxmin + d))
        xmax = int(max(subxmax, subxmax + d))

        # Gaussians in large image
        npix = np.prod(subshape)
        npixperpeak = 2.0
        npeaks = int(npix / npixperpeak)
        sx = npixperpeak * 1.5
        margin = int(10 * sx)
        xmin -= margin
        xmax += margin

        shape = subshape
        if transposed:
            shape = (xmax - xmin + 1, 1)
        else:
            shape = (1, xmax - xmin + 1)

        np.random.seed(1)

        if realistic:
            x0, npeaks = genpos1d(xmin, xmax, npeaks)
        else:
            x0, npeaks = genpos1d(subxmin, subxmax, npeaks)
        ind = (x0 > 5) | (x0 < 5)
        x0 = x0[ind]
        npeaks = sum(ind)
        sx = random(sx * 0.8, sx * 1.2, npeaks)
        A = random(1, 5, npeaks)
        data = tuple(zip(x0, sx, A))

        # Plot full images
        if plot:
            xv = np.arange(xmin, xmax + 1)
            img = gettransformedvector(xv, data).reshape(shape)
            import matplotlib.pyplot as plt

            plt.figure(2)
            plt.plot(xv.flat, img.flat)
            plt.show()

        # Stack of transformed subimages
        if realistic:
            s = subshape
            xv = np.arange(subxmin, subxmax + 1)
        else:
            s = shape
            xv = np.arange(xmin, xmax + 1)

        ret = np.empty((nimages,) + s, dtype=np.float32)
        xy = np.copy(xv)
        ret[0] = gettransformedvector(xy, data).reshape(s)
        for i in range(1, nimages):
            xy = xy + Mcof[1 - dimi, 2]
            ret[i] = gettransformedvector(xy, data).reshape(s)
    else:
        # Shape of a subimage
        subshape = (ndim1, ndim2)
        subxmin = -subshape[1] // 2
        subymin = -subshape[0] // 2
        subxmax = subshape[1] // 2
        subymax = subshape[0] // 2
        subxmax = subshape[1] - 1
        subymax = subshape[0] - 1
        subshape = (subymax - subymin + 1, subxmax - subxmin + 1)

        # Determine shape of large image (transform corners)
        xy = np.empty((3, 4))
        xy[0, :] = [subxmin, subxmax, subxmin, subxmax]
        xy[1, :] = [subymin, subymin, subymax, subymax]
        xy[2, :] = [1, 1, 1, 1]
        myminmax = np.append(np.min(xy, axis=1), np.max(xy, axis=1))

        for i in range(1, nimages):
            xy = np.dot(Mcoord, xy)
            xy[0, :] /= xy[2, :]
            xy[1, :] /= xy[2, :]
            myminmax[0:3] = np.minimum(myminmax[0:3], np.min(xy, axis=1))
            myminmax[3:] = np.maximum(myminmax[3:], np.max(xy, axis=1))
        xmin = int(myminmax[0])
        ymin = int(myminmax[1])
        xmax = int(myminmax[3])
        ymax = int(myminmax[4])

        # Gaussians in large image
        npix = np.prod(subshape)
        npixperpeak = 10.0
        npeaks = int(npix / npixperpeak)
        sxy = np.sqrt(npixperpeak)
        margin = int(10 * sxy)
        xmin -= margin
        ymin -= margin
        xmax += margin
        ymax += margin
        shape = (ymax - ymin + 1, xmax - xmin + 1)

        np.random.seed(1)
        if realistic:
            x0, y0, npeaks = genpos2d(xmin, xmax, ymin, ymax, npeaks)
        else:
            x0, y0, npeaks = genpos2d(subxmin, subxmax, subymin, subymax, npeaks)
        sx = random(sxy * 0.8, sxy * 1.2, npeaks)
        sy = random(sxy * 0.8, sxy * 1.2, npeaks)
        rho = random(0, 0.2, npeaks)
        A = random(1, 5, npeaks)
        data = tuple(zip(x0, y0, sx, sy, rho, A))

        # Plot full image
        if plot:
            xv, yv = np.meshgrid(np.arange(xmin, xmax + 1), np.arange(ymin, ymax + 1))
            xv = xv.reshape((1, shape[0] * shape[1]))
            yv = yv.reshape((1, shape[0] * shape[1]))
            xy = np.vstack((xv, yv, np.ones_like(xv)))
            img = gettransformedimage(xv, yv, data).reshape(shape)
            import matplotlib.pyplot as plt

            plt.figure(2)
            plt.subplot(111)
            plt.imshow(img, origin="lower", interpolation="nearest")
            plt.show()

        # Stack of transformed subimages
        if realistic:
            s = subshape
            xv, yv = np.meshgrid(
                np.arange(subxmin, subxmax + 1), np.arange(subymin, subymax + 1)
            )
            ox, oy = subxmin, subymin
        else:
            s = shape
            xv, yv = np.meshgrid(np.arange(xmin, xmax + 1), np.arange(ymin, ymax + 1))
            ox, oy = xmin, ymin

        ret = np.empty((nimages,) + s, dtype=np.float32)
        xv = xv.reshape((1, s[0] * s[1]))
        yv = yv.reshape((1, s[0] * s[1]))
        xy = np.vstack((xv, yv, np.ones_like(xv)))
        ret[0] = gettransformedimage(xv, yv, data).reshape(s)
        for i in range(1, nimages):
            xy = np.dot(Mcof, xy)  # coordinates from new frame to old frame
            xy[0, :] /= xy[2, :]
            xy[1, :] /= xy[2, :]
            ret[i] = gettransformedimage(xy[0, :], xy[1, :], data).reshape(s)

        # Relative change-of-frame in subimage pixel frame
        C = np.identity(3, dtype=Mcof.dtype)
        Cinv = np.identity(3, dtype=Mcof.dtype)
        C[0:2, 2] = [ox, oy]  # image pixel frame to subimage pixel frame
        Cinv[0:2, 2] = -C[0:2, 2]
        Mcof = np.dot(np.dot(Cinv, Mcof), C)
        Mcoord = np.dot(np.dot(Cinv, Mcoord), C)

    # Mcoord: change-of-frame matrix of the back-transformation
    if inverse:
        ret = ret.max() - ret
    return [ret.copy() for _ in range(nstacks)], Mcoord, stackdim
