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


from ..utils import listtools
from ..utils import instance
from ..utils import units

import numpy as np
from scipy.interpolate import interp1d


class discrete(object):

    def __init__(self, lines, intensities=None):
        """
        Args:
            lines(ureg.Quantity): in keV(default), nm, ...
            intensities(array): line intensities
        """

        lines = units.Quantity(lines, "keV")
        unit = lines.units
        lines = lines.magnitude

        self.bnumber = not instance.isarray(lines)

        if instance.isarray(lines):
            self.bnumber = False
            self.size = len(lines)
        else:
            self.bnumber = True
            self.size = 1

        # Intensities
        if intensities is None:
            if self.bnumber:
                intensities = 1.
            else:
                intensities = [1./self.size]*self.size
        else:
            if not instance.isarray(intensities) and not self.bnumber:
                intensities = [intensities]

        # Handle duplicates
        if self.bnumber:
            self.intensities = intensities
        else:
            # Sum equal lines
            lines, intensities = listtools.sumrepeats(lines, intensities)

            # Sort by line
            lines, self.intensities = listtools.sort2lists(lines, intensities)

        # Lines
        self.lines = units.Quantity(lines, unit)

    def __getstate__(self):
        return {'lines': self.lines,
                'intensities': self.intensities,
                'bnumber': self.bnumber,
                'size': self.size}

    def __setstate__(self, state):
        self.lines = state['lines']
        self.intensities = state['intensities']
        self.bnumber = state['bnumber']
        self.size = state['size']

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return np.array_equal(self.lines, other.lines) and \
                np.array_equal(self.intensities, other.intensities) and \
                self.bnumber == other.bnumber and \
                self.size == other.size
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __add__(self, val):
        if isinstance(val, self.__class__):
            lines1 = self.lines.magnitude
            intensities1 = self.intensities
            if self.bnumber:
                lines1 = [lines1]
                intensities1 = [intensities1]

            lines2 = val.lines.to(self.lines.units).magnitude
            intensities2 = val.intensities
            if val.bnumber:
                lines2 = [lines2]
                intensities2 = [intensities2]

            if isinstance(lines1, np.ndarray):
                lines1 = lines1.tolist()
            if isinstance(lines2, np.ndarray):
                lines2 = lines2.tolist()

            return self.__class__(units.Quantity(lines1 + lines2, self.lines.units),
                                  intensities1 + intensities2)
        else:
            raise NotImplementedError

    @property
    def total(self):
        if self.bnumber:
            return self.intensities
        else:
            return sum(self.intensities)

    @property
    def ratios(self):
        if self.bnumber:
            return 1
        else:
            return np.asarray(self.intensities)/float(self.total)

    @property
    def energies(self):
        return self.lines.to('keV', 'spectroscopy').magnitude

    @property
    def nlines(self):
        if self.bnumber:
            return 1
        else:
            return len(self.intensities)

    def weightedsum(self, x):
        if self.bnumber:
            return x
        else:
            x = instance.asarray(x)
            return sum(x*self.ratios)

    def sample(self, x, y):
        x = units.Quantity(x, 'keV').to('keV', 'spectroscopy').magnitude
        x, y = listtools.sort2lists(x, y)
        f = interp1d(x, y, kind='linear', bounds_error=False,
                     fill_value=(y[0], y[-1]))
        return self.weightedsum(f(self.energies))


class dark(object):

    def sample(self, x, y):
        return 0
