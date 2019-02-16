# -*- coding: utf-8 -*-
#
#   Copyright (C) 2017 European Synchrotron Radiation Facility, Grenoble, France
#
#   Principal author:   Wout De Nolf (wout.de_nolf@esrf.eu)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

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
