import time
import os
import numpy as np
from .utils import mkdir
from .spec import spec


def array2SpecMca(data):
    """Write a python array into a Spec array.
    Return the string containing the Spec array
    """
    tmpstr = "@A "
    length = len(data)
    for idx in range(0, length, 16):
        if idx + 15 < length:
            for i in range(0, 16):
                tmpstr += "%.8g " % data[idx + i]
            if idx + 16 != length:
                tmpstr += "\\"
        else:
            for i in range(idx, length):
                tmpstr += "%.8g " % data[i]
        tmpstr += "\n"
    return tmpstr


def save(mca, filename, mode="w", zero=0, gain=1, ch0=None, ch1=None):
    # http://www.esrf.eu/blissdb/macros/macdoc.py?macname=saveload.mac
    path = os.path.dirname(filename)
    if path:
        mkdir(path)
    if ch0 is None:
        ch0 = 0
    if ch1 is None:
        ch1 = ch0 + len(mca) - 1
    with open(filename, mode) as ffile:
        ffile.write("#F %s\n" % filename)
        ffile.write("#D %s\n" % (time.ctime(time.time())))
        ffile.write("\n")
        ffile.write("#S 1 %s\n" % "spectrocrunch")
        ffile.write("#D %s\n" % (time.ctime(time.time())))
        ffile.write("#@MCA %16C\n")
        # 1: data reduction coefficient???
        ffile.write("#@CHANN %d %d %d 1\n" % (ch1 - ch0 + 1, ch0, ch1))
        ffile.write("#@CALIB %.7g %.7g %.7g\n" % (zero, gain, 0))
        # ffile.write("#@CTIME 1 0.9 0.99")  # preset, live, real
        ffile.write(array2SpecMca(mca))
        ffile.write("\n")


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
    ch0 = ch1 = 0
    for s in h:
        if s.startswith("#@CALIB"):
            coeff = s.replace("#@CALIB", "").strip()
            coeff = np.array(map(float, coeff.split(" ")))
        elif s.startswith("#@CHANN"):
            tmp = s.replace("#@CHANN", "").strip()
            _, ch0, ch1, _ = list(map(int, tmp.split(" ")))
    if len(mca):
        # Bug in MCA reading: last channel is not read
        mca = np.append(mca, 0)
    if ch0 == ch1:
        ch0 = 0
        ch1 = len(mca) - 1
    channels = np.arange(ch0, ch1 + 1)
    energy = sum(c * channels**i for i, c in enumerate(coeff))
    return mca, channels, energy, coeff
