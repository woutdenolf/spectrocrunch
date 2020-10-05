# -*- coding: utf-8 -*-

import os, sys

sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import scipy.ndimage.interpolation


def pad(img, extend):
    """Apply padding"""
    pad = (
        (max(extend[0][0], 0), max(extend[0][1], 0)),
        (max(extend[1][0], 0), max(extend[1][1], 0)),
    )
    if np.count_nonzero(pad) != 0:
        return np.pad(img, pad, "constant", constant_values=-2)
    else:
        return img


def crop(img, extend):
    """Apply cropping"""
    dim1, dim2 = img.shape
    crop = (
        (max(-extend[0][0], 0), dim1 - max(-extend[0][1], 0)),
        (max(-extend[1][0], 0), dim2 - max(-extend[1][1], 0)),
    )
    if crop[0][0] != 0 or crop[1][0] != 0 or crop[0][1] != dim1 or crop[1][1] != dim2:
        return img[crop[0][0] : crop[0][1], crop[1][0] : crop[1][1]]
    else:
        return img


def extendfromtransformation(offsets, shape, pad=True):
    """If all transformations are known, padding/cropping can be calculated"""

    extend = ((0, 0), (0, 0))
    n1, n2 = shape

    o1min = np.floor(np.min(offsets[:, 0])).astype(np.int)
    o2min = np.floor(np.min(offsets[:, 1])).astype(np.int)
    o1max = np.ceil(np.max(offsets[:, 0])).astype(np.int)
    o2max = np.ceil(np.max(offsets[:, 1])).astype(np.int)

    if pad:
        extend = ((o1max, -o1min), (o2max, -o2min))
    else:
        extend = ((o1min, -o1max), (o2min, -o2max))

    return extend


def transform(img, extend, offset):
    img = pad(img, extend)
    img = scipy.ndimage.interpolation.affine_transform(
        img, np.identity(2), offset=offset, cval=-1
    )
    img = crop(img, extend)
    return img


if __name__ == "__main__":

    mi = -1.1
    ma = 2
    bpad = False

    offsets = np.zeros((2, 2), dtype=np.float)
    offsets[0, 0] = mi
    offsets[0, 1] = mi
    offsets[1, 0] = ma
    offsets[1, 1] = ma

    img = np.arange(5 * 6).reshape((5, 6))

    extend = extendfromtransformation(offsets, img.shape, pad=bpad)
    print(extend)

    imga = transform(img, extend, (mi, mi))
    imgb = transform(img, extend, (ma, ma))
    print(imga)
    print(imgb)
    if bpad:
        tmp1 = imga < 0
        tmp2 = imgb < 0
        imga[tmp1] = 0
        imgb[tmp1] = 0
        imga[tmp2] = 0
        imgb[tmp2] = 0
        print(imga)
        print(imgb)
