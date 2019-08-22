# -*- coding: utf-8 -*-

from ..utils.classfactory import with_metaclass
from . import base
import numpy as np


class DiodeGeometry(with_metaclass(base.SolidAngle)):
    pass


class SXM_IODET(DiodeGeometry):
    aliases = ['iodet', 'iodet1', 'iodet2', 'sxm_iodet1', 'sxm_iodet2']

    def __init__(self, **kwargs):
        kwargs["anglein"] = kwargs.pop("anglein", 90)
        kwargs["angleout"] = kwargs.pop("angleout", 70)
        kwargs["azimuth"] = kwargs.pop("azimuth", 0)
        kwargs["solidangle"] = kwargs.pop("solidangle", 4*np.pi*0.4)
        super(SXM_IODET, self).__init__(**kwargs)


class SXM_FDET(DiodeGeometry):
    aliases = ['fdet']

    def __init__(self, **kwargs):
        kwargs["anglein"] = kwargs.get("anglein", 62)
        kwargs["angleout"] = kwargs.get("angleout", 58)  # 49
        kwargs["azimuth"] = kwargs.get("azimuth", 0)
        kwargs["solidangle"] = kwargs.pop("solidangle", 4*np.pi*0.4)
        super(SXM_FDET, self).__init__(**kwargs)


class SXM_IDET(DiodeGeometry):
    aliases = ['idet']

    def __init__(self, **kwargs):
        kwargs["anglein"] = kwargs.get("anglein", 90)
        kwargs["angleout"] = kwargs.get("angleout", -90)
        kwargs["azimuth"] = kwargs.get("azimuth", 0)
        kwargs["solidangle"] = kwargs.pop("solidangle", 4*np.pi*0.4)
        super(SXM_IDET, self).__init__(**kwargs)


class ID16B_IC(DiodeGeometry):
    aliases = ['id16b_IC']

    def __init__(self, **kwargs):
        kwargs["anglein"] = kwargs.get("anglein", 90)
        kwargs["angleout"] = kwargs.get("angleout", -90)
        kwargs["azimuth"] = kwargs.get("azimuth", 0)
        kwargs["solidangle"] = kwargs.pop("solidangle", 4*np.pi*0.4)
        super(ID16B_IC, self).__init__(**kwargs)


factory = DiodeGeometry.factory
registry = DiodeGeometry.clsregistry
