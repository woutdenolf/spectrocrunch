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
from ..utils import instance
from ..utils import units
from ..utils import lut

import numpy as np
import logging

logger = logging.getLogger(__name__)


class Discrete(lut.LUT):

    def __init__(self, lines=None, intensities=None, **kwargs):
        """
        Args:
            lines(Quantity): in keV(default), nm, ...
            intensities(num or array): line intensities
        """
        if lines is None:
            kwargs.pop('x', None)
        else:
            lines = units.Quantity(lines, "keV")
            kwargs['x'] = lines
        if intensities is None:
            kwargs.pop('y', None)
        else:
            intensities = units.Quantity(intensities, "dimensionless")
            kwargs['y'] = intensities
        kwargs['default'] = units.Quantity(0, "dimensionless")
        super(Discrete, self).__init__(**kwargs)

    @property
    def lines(self):
        return self.x

    @property
    def intensities(self):
        return self.y.magnitude

    @property
    def total(self):
        return sum(self.intensities)

    @property
    def ratios(self):
        return self.intensities/float(self.total)

    @property
    def energies(self):
        return self.lines.to('keV', 'spectroscopy').magnitude

    @property
    def nlines(self):
        return len(self)

    def sample(self, lut):
        if self.isempty():
            y = lut.default
        else:
            lines = self.lines
            if not lut.isempty():
                lines = lines.to(lut.xunits, 'spectroscopy')
            y = lut(lines)
            y = sum(y*self.ratios)
        return y


class Dark(Discrete):

    def __init__(self, **kwargs):
        super(Dark, self).__init__(lines=None, intensities=None, **kwargs)

    def __add__(self, xy):
        logger.warning('No lines can be added to a Dark spectrum')

    def __iadd__(self, xy):
        logger.warning('No lines can be added to a Dark spectrum')
