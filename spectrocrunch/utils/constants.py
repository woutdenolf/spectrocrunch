# -*- coding: utf-8 -*-

from ..patch.pint import ureg
from ..utils import units


def temperatureinkelvin(T):
    """
    Args:
        T(num|array): temperature in deg C
    Returns:
        num|array: keV
    """
    T = units.Quantity(T, ureg.degC)
    return T.to(ureg.kelvin)


def eholepair_si(T=21):
    """
    Args:
        T(num|array): temperature in deg C
    Returns:
        num|array: keV
    """
    # https://doi.org/10.1016/j.nima.2007.03.020
    T = temperatureinkelvin(T)
    x = units.Quantity([80, 270], ureg.kelvin)  # K
    y = units.Quantity([3.77, 3.68], "eV")

    m = (y[1] - y[0]) / (x[1] - x[0])
    b = y[1] - m * x[1]

    return (m * T + b).to("eV")
