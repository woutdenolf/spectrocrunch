"""
Statistical methods.
"""

import numpy as np


def outlierdetection(x, nsigma, pdf="normal", noutliers=0):

    # Detect outliers using the median-absolute-deviation
    # MAD = cte*median(x-median(x))
    # |(x-median(x))/MAD|> nsigma

    if x.size == 0:
        np.full(1, False, dtype=bool)

    # Deviation form medium
    diff = abs(x - np.median(x))

    # Fixed number of outliers
    if noutliers != 0:
        ret = np.full(x.size, False, dtype=bool)
        ret[(-diff).argsort()[0:noutliers]] = True
        return ret

    # median-absolute-deviation
    if pdf == "normal":
        MAD = 1.4826 * np.median(diff)
    else:
        MAD = np.median(diff)

    # outliers
    if MAD == 0:
        return np.full(x.shape, False, dtype=bool)
    else:
        return diff / MAD > nsigma
