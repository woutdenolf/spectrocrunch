# -*- coding: utf-8 -*-

import numpy as np
from .spec import spec


def read(filename):
    """
    Args:
        filename(str)
    Returns:
        tuple(array): mca, channels, energy, energy calibration coefficients
    """
    f = spec(filename)
    mca = f.getdata2(1, [])
    h = f.getKeyInfo("1.1")["Header"]
    coeff = np.arange(2)
    try:
        i = h.index("$MCA_CAL:")
    except ValueError:
        pass
    else:
        coeff = np.array(list(map(float, h[i + 2].split(" "))))
    channels = np.arange(len(mca))
    energy = sum(c * channels ** i for i, c in enumerate(coeff))
    return mca, channels, energy, coeff
