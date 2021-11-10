# -*- coding: utf-8 -*-

import sys
import difflib
import numpy as np
import pandas as pd
from ..utils import units
from ..utils import instance
from ..utils import listtools


def arg_closest_num(arr, value, check):
    if check(value):
        return value
    else:
        return np.argmin(np.abs(arr - value))


def arg_closest_string(arr, value, check):
    if check(value):
        return value
    else:
        s = difflib.SequenceMatcher()
        s.set_seq2(value)
        j, ratio = 0, 0
        for i, x in enumerate(arr):
            s.set_seq1(x)
            iratio = s.ratio()
            if iratio > ratio:
                j, ratio = i, iratio
        return j


def arg_index(arr, value, check):
    if check(value):
        return value
    else:
        return arr.index(value)


class Axis(object):
    def __init__(self, params, type=None, name=None, title=None, precision=None):
        if type is None:
            type = "quantitative"
            # quantitative: number
            # nominal: unordered categorical
            # ordinal: ordered categorical
            # temporal: date/time
        self.type = type
        self.name = name
        self.title = title
        self.values = params
        self.precision = precision

    def __getitem__(self, index):
        return self.values[index]

    def __setitem__(self, index, values):
        if self.type == "quantitative":
            values = units.Quantity(values, units=self.units)
            self.values[index] = values
        else:
            self.values.iloc[index] = values

    @property
    def initargs(self):
        if isinstance(self._params, tuple):
            return self._params
        else:
            return (self._params,)

    @property
    def initkwargs(self):
        return {
            "type": self.type,
            "name": self.name,
            "title": self.title,
            "precision": self.precision,
        }

    @property
    def units(self):
        if self.type == "quantitative":
            return self.values.units
        else:
            return None

    @units.setter
    def units(self, value):
        if self.type != "quantitative":
            raise RuntimeError("{} axis has no units".format(self.type))
        self.values.ito(value)

    @property
    def unit_suffix(self):
        if self.units:
            return "({:~})".format(self.units)
        else:
            return ""

    @property
    def name(self):
        if self._name:
            return self._name
        else:
            return self.__class__.__name__

    @name.setter
    def name(self, value):
        self._name = value

    @property
    def title_nounits(self):
        if self._title:
            title = self._title
        else:
            title = self.name
        suffix = self.unit_suffix
        if suffix:
            title = title.replace(suffix, "")
        return title

    @property
    def title(self):
        if self._title:
            title = self._title
        else:
            title = self.name
        suffix = self.unit_suffix
        if suffix:
            if suffix not in title:
                title = "{} {}".format(title, suffix)
        return title

    @title.setter
    def title(self, value):
        self._title = value

    @property
    def magnitude(self):
        """np.ndarray or pandas.Index"""
        if self.type == "quantitative":
            return self.values.magnitude
        else:
            return self.values

    def umagnitude(self, u):
        if self.type == "quantitative":
            return self.values.to(u).magnitude
        else:
            return self.values

    @property
    def values(self):
        return self._values

    @values.setter
    def values(self, params):
        """pint.Quantity or pandas.Index"""
        self._params = params
        self._values = params
        if self.type == "quantitative":
            self._values = units.asqarray(self._values)

    @property
    def dtype(self):
        return self.values.dtype

    @property
    def start(self):
        return self.values[0]

    @property
    def end(self):
        return self.values[-1]

    @property
    def size(self):
        return len(self.values)

    @property
    def nsteps(self):
        return self.size - 1

    @property
    def precision(self):
        p = self._precision
        if self.type == "quantitative":
            p = p.to(self.units)
        return p

    @precision.setter
    def precision(self, value):
        if value is None:
            value = 0
        if self.type == "quantitative":
            value = units.Quantity(value, units=self.units)
        self._precision = value

    def __len__(self):
        return self.size

    def __str__(self):
        return self.title

    def __repr__(self):
        return "{}({})".format(self.name, len(self))

    def __eq__(self, other):
        if self.size != other.size:
            return False
        if self.type != "quantitative":
            return all(self.values == other.values)
        # threshold = max(self.precision, other.precision.to(self.units)).magnitude
        threshold = self.precision.magnitude
        values = self.magnitude
        uvalues = other.umagnitude(self.units)

        mask_nan = np.isnan(values)
        umask_nan = np.isnan(uvalues)
        if not (mask_nan == umask_nan).all():
            return False
        if mask_nan.all():
            return True

        mask_not_nan = ~mask_nan
        values = values[mask_not_nan]
        uvalues = uvalues[mask_not_nan]
        diff = max(abs(values - uvalues))
        return diff <= threshold

    def __ne__(self, other):
        return not self.__eq__(other)

    def _extract_magnitude(self, values):
        if instance.isquantity(values):
            return units.umagnitude(values, units=self.units)
        elif isinstance(values, Axis):
            return values.umagnitude(self.units)
        else:
            return values

    def locate(self, values, detectindex=False):
        """
        Get indices of values:
         None: entire range slice(None)
         string: position of closest string
         integer and detectindex==True: position of closest value when Axis are not integers
         float: position closest value

        Args:
            values(array or num or None):
            detectindex(bool):
        Returns:
            list(int)|int: indices of values
        """
        if values is None:
            return slice(None)

        # Axis values (always an array)
        xold = self.magnitude
        if self.type == "quantitative":
            func = arg_closest_num
        elif instance.isstring(xold[0]):
            func = arg_closest_string
        else:
            func = arg_index
            try:
                xold = xold.tolist()
            except AttributeError:
                pass

        def isindex(v):
            return False

        if detectindex:
            if self.type == "quantitative":
                if instance.isinteger(xold[0]):

                    def isindex(v):
                        return instance.isinteger(v)

        xnew = self._extract_magnitude(values)
        if instance.isarray(xnew):
            if instance.isstring(xnew[0]):
                func = arg_closest_string
                xold = list(map(str, xold))
            elif instance.isstring(xold[0]):
                xnew = list(map(str, xnew))
            return [func(xold, v, isindex) for v in xnew]
        else:
            if instance.isstring(xnew):
                func = arg_closest_string
                xold = list(map(str, xold))
            return func(xold, xnew, isindex)

    def interpolate(self, values):
        """Get old and new axes values for interpolation"""
        xold = self.magnitude
        xnew = self._extract_magnitude(values)
        if self.type == "quantitative":
            if xold[1] < xold[0]:
                # Reverse axis
                m = np.nanmax(xold)
                xold = m - xold
                if xnew is not None:
                    xnew = m - xnew
            if xnew is None:
                return xold, xold
            else:
                return xold, xnew
        else:
            # Use index for interpolation
            ind = range(len(self))
            if xnew is None:
                return ind, ind
            else:
                indnew = []
                if not instance.isarray(xnew):
                    xnew = [xnew]
                if instance.isarray(xold):
                    try:
                        xold = xold.tolist()
                    except AttributeError:
                        pass
                else:
                    xold = [xold]
                for x in xnew:
                    indnew.append(xold.index(x))
                return ind, indnew

    def to(self, u):
        if self.type != "quantitative":
            raise RuntimeError("{} axis has no units".format(self.type))
        ret = self.__class__(*self.initargs, **self.initkwargs)
        ret.units = u
        return ret

    def simplify(self):
        ax = self
        if self.type == "quantitative":
            kwargs = {}
            if self.size == 1:
                return AxisNumber(self.start, **self.initkwargs)
            elif self.size == 2:
                return AxisRegular(self.start, self.end, 1, **self.initkwargs)
            else:
                diff = abs(np.diff(np.diff(self.magnitude)))
                if max(diff) <= self.precision.magnitude:
                    return AxisRegular(
                        self.start, self.end, self.size - 1, **self.initkwargs
                    )
                else:
                    n = len(self)
                    ind = np.array(
                        [0]
                        + (np.where(diff > self.precision.magnitude)[0] + 1).tolist()
                        + [n - 1]
                    )
                    nseg = len(ind) - 1
                    if nseg < n / 3.0:
                        # Limit the number of segments to 1/3 of the number of values
                        limits = self.values[ind]
                        nsteps = np.diff(ind)
                        return AxisSegments(limits, nsteps, **self.initkwargs)

        return self

    @property
    def start_stepsize(self):
        return self[1] - self[0]

    @property
    def end_stepsize(self):
        return self[-1] - self[-2]

    @classmethod
    def factory(cls, data, **kwargs):
        if not isinstance(data, Axis):
            data = Axis(data, **kwargs)
        return data.simplify()

    def newstart(self, newstart):
        """Expand/crop from the start (make sure newstart is included)

        Args:
            newstart(num)

        Returns:
            Axis|None
        """
        x = self.magnitude
        stepsize = self.start_stepsize
        newstart = units.Quantity(newstart, units=self.units)
        oldstart = self.start
        nstepsadd = (oldstart - newstart) / stepsize
        nstepsadd = int(np.ceil(nstepsadd.to("dimensionless").magnitude))
        if nstepsadd > 0:
            add = oldstart - np.arange(1, nstepsadd + 1)[::-1] * stepsize
            x = np.append(add.magnitude, x)
        elif nstepsadd < 0:
            x = x[-nstepsadd:]
        if len(x):
            return self.factory(units.Quantity(x, units=self.units))
        else:
            return None

    def newend(self, newend):
        """Expand/crop from the start (make sure newend is included)

        Args:
            newend(num)

        Returns:
            Axis|None
        """
        x = self.magnitude
        stepsize = self.end_stepsize
        newend = units.Quantity(newend, units=self.units)
        oldend = self.end
        nstepsadd = (newend - oldend) / stepsize
        nstepsadd = int(np.ceil(nstepsadd.to("dimensionless").magnitude))
        if nstepsadd > 0:
            add = oldend + np.arange(1, nstepsadd + 1) * stepsize
            x = np.append(x, add.magnitude)
        elif nstepsadd < 0:
            x = x[:nstepsadd]
        if len(x):
            return self.factory(units.Quantity(x, units=self.units))
        else:
            return None

    def newlimits(self, newstart, newend):
        """Expand/crop from the start (make sure limits are included)

        Args:
            newstart(num)
            newend(num)

        Returns:
            Axis|None
        """
        axis = self.newstart(newstart)
        if axis is not None:
            axis = axis.newend(newend)
        return axis


factory = Axis.factory


class _AxisRegular(Axis):
    def __init__(self, *params, **kwargs):
        kwargs["type"] = "quantitative"
        super(_AxisRegular, self).__init__(params, **kwargs)

    def __setitem__(self, ind, values):
        raise RuntimeError("Axis values are not editable")

    @property
    def start(self):
        return self._start

    @start.setter
    def start(self, value):
        self._start = value

    @property
    def end(self):
        return self._end

    @end.setter
    def end(self, value):
        self._end = value

    @property
    def stepsize(self):
        return self._stepsize

    @stepsize.setter
    def stepsize(self, value):
        self._stepsize = value

    @property
    def start_stepsize(self):
        return self.stepsize

    @property
    def end_stepsize(self):
        return self.stepsize

    @property
    def size(self):
        return self._size

    @size.setter
    def size(self, value):
        self._size = value

    @property
    def nsteps(self):
        return self._size - 1

    @nsteps.setter
    def nsteps(self, value):
        self._size = value + 1

    @property
    def units(self):
        return self.start.units

    @units.setter
    def units(self, value):
        self.values.ito(value)
        self.start.ito(value)
        self.end.ito(value)
        self.stepsize.ito(value)

    def __repr__(self):
        return "{}(start={:~},end={:~},step={:~},size={})".format(
            self.name, self.start, self.end, self.stepsize, len(self)
        )


class AxisRegular(_AxisRegular):
    """start, end, nsteps"""

    @_AxisRegular.start.setter
    def start(self, value):
        self.values = value, self.end, self.size

    @_AxisRegular.end.setter
    def end(self, value):
        self.values = self.start, value, self.size

    @_AxisRegular.size.setter
    def size(self, value):
        self.values = self.start, self.end, value - 1

    @_AxisRegular.nsteps.setter
    def nsteps(self, value):
        self.values = self.start, self.end, value

    @_AxisRegular.values.setter
    def values(self, params):
        self._params = params
        self._start = units.Quantity(params[0])
        u = self._start.units
        self._end = units.Quantity(params[1], units=u).to(u)
        self._size = params[2] + 1
        if self.size == 1:
            self._stepsize = units.Quantity(0, units=u).to(u)
        else:
            self._stepsize = (self.end - self.start) / self.nsteps
        self._values = units.Quantity(
            np.linspace(self.start.magnitude, self.end.magnitude, self.size), units=u
        )


class AxisRegularInc(_AxisRegular):
    """start, stepsize, nsteps"""

    @_AxisRegular.start.setter
    def start(self, value):
        self.values = value, self.stepsize, self.size

    @_AxisRegular.stepsize.setter
    def stepsize(self, value):
        self.values = self.start, value, self.size

    @_AxisRegular.size.setter
    def size(self, value):
        self.values = self.start, self.stepsize, value - 1

    @_AxisRegular.nsteps.setter
    def nsteps(self, value):
        self.values = self.start, self.stepsize, value

    @_AxisRegular.values.setter
    def values(self, params):
        self._params = params
        self._start = units.Quantity(params[0])
        u = self._start.units
        self._stepsize = units.Quantity(params[1], units=u).to(u)
        self._size = params[2] + 1
        self._end = self.start + self.stepsize * self.nsteps
        self._values = np.arange(self.size) * self.stepsize + self.start


def arange(*args):
    return factory(np.arange(*args))


def zapscan(start, end, npixels, unit=None, **kwargs):
    pixelsize = (end - start) / float(npixels)
    if unit is None:
        unit = "dimensionless"
    start = units.Quantity(start + pixelsize / 2, units=unit)
    return AxisRegularInc(start, pixelsize, npixels - 1, **kwargs)


def ascan(start, end, intervals, unit=None, **kwargs):
    if unit is None:
        unit = "dimensionless"
    start = units.Quantity(start, units=unit)
    return AxisRegular(start, end, intervals, **kwargs)


class AxisNumber(_AxisRegular):
    @_AxisRegular.start.setter
    def start(self, value):
        self.values = (value,)

    @_AxisRegular.values.setter
    def values(self, params):
        self._params = params
        self._start = units.Quantity(params[0])
        u = self._start.units
        self._end = self.start
        self._size = 1
        self._stepsize = units.Quantity(0, units=u)
        self._values = units.Quantity([self.start.magnitude], units=u)

    def __repr__(self):
        return "{}(={:~})".format(self.name, self.start)


class AxisSegments(Axis):
    def __init__(self, *params, **kwargs):
        kwargs["type"] = "quantitative"
        super(AxisSegments, self).__init__(params, **kwargs)

    def __setitem__(self, ind, values):
        raise RuntimeError("Axis values are not editable")

    @property
    def start(self):
        return self.limits[0]

    @property
    def end(self):
        return self.limits[-1]

    @property
    def limits(self):
        return self._limits

    @limits.setter
    def limits(self, values):
        self.values = values, self.nsteps

    @property
    def stepsizes(self):
        return self._stepsizes

    @stepsizes.setter
    def stepsizes(self, values):
        limits = self.limits

        u = units.Quantity(values[0]).units
        if u == "dimensionless":
            u = limits.units

        def func(x):
            return units.Quantity(x, units=u)

        nsteps = []
        for i, stepsize in enumerate(values):
            stepsize = func(stepsize)
            a, b = limits[i], limits[i + 1]
            nsteps.append(int(round((b - a) / stepsize)))

        self.values = limits, nsteps

    @property
    def start_stepsize(self):
        return self.stepsizes[0]

    @property
    def end_stepsize(self):
        return self.stepsizes[1]

    @property
    def nsteps(self):
        return self._nsteps

    @nsteps.setter
    def nsteps(self, values):
        self.values = self.limits, values

    @staticmethod
    def mergeargs(limits, nsteps):
        nsteps = np.asarray(nsteps).tolist() + [None]
        return tuple(list(listtools.flatten(zip(limits, nsteps)))[:-1])

    @property
    def nsegments(self):
        return len(self._nsteps)

    @Axis.values.setter
    def values(self, params):
        limits, nsteps = params
        if len(limits) != len(nsteps) + 1:
            raise ValueError("Number of segments does not match the number of limits")
        self._params = params

        u = units.Quantity(limits[0]).units

        def func(x):
            return units.Quantity(x, units=u).to(u).magnitude

        lmts = []
        stepsizes = []
        values = []
        for i in range(len(nsteps)):
            start, end = func(limits[i]), func(limits[i + 1])
            nspt = func(nsteps[i])
            inc = (end - start) / float(nspt)
            values += (start + inc * np.arange(nspt)).tolist()
            lmts.append(start)
            stepsizes.append(inc)
        lmts.append(end)
        values.append(end)

        self._stepsizes = units.Quantity(stepsizes, units=u)
        self._limits = units.Quantity(lmts, units=u)
        self._nsteps = np.asarray(nsteps)
        self._values = units.Quantity(values, units=u)

    def __repr__(self):
        s = "".join(
            list(
                "{:~}--({}x{:~})--".format(lim, n, step)
                for lim, step, n in zip(self.limits, self.stepsizes, self.nsteps)
            )
        )
        s = "{}{:~}".format(s, self.limits[-1])
        return "{}({})".format(self.name, s)
