# -*- coding: utf-8 -*-

from ..utils.classfactory import with_metaclass
from . import base

import numpy as np


class Geometry(with_metaclass(base.FlatSample)):

    def __init__(self, **kwargs):
        super(Geometry, self).__init__(**kwargs)


class Perpendicular(Geometry):

    def __init__(self, **kwargs):
        super(Perpendicular, self).__init__(anglein=90,
                                            angleout=-90,
                                            **kwargs)


factory = Geometry.factory
registry = Geometry.clsregistry
