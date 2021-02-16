# -*- coding: utf-8 -*-

from contextlib import contextmanager
import numpy as np

try:
    from contextlib import ExitStack
except ImportError:
    from contextlib2 import ExitStack

from .regulargrid import RegularGrid
from . import axis
from . import nxresult
from ..utils import indexing
from ..io import nxfs


class NXSignalRegularGrid(RegularGrid):
    def __init__(self, signal, stackdim=None):
        nxdata = signal.parent
        axes = [
            axis.factory(values, name=name, title=attrs["title"], type="quantitative")
            for name, values, attrs in nxdata.axes
        ]
        self.signal = signal
        super(NXSignalRegularGrid, self).__init__(axes, None, stackdim=stackdim)

        stackdim_prev = nxdata.stackdim(default=self.DEFAULT_STACKDIM)
        stackdim = self.stackdim
        if stackdim <= stackdim_prev:
            stackdim_prev += 1
        self._stack_dims = [self.stackdim, stackdim_prev]

    @property
    def stack_dims(self):
        return self._stack_dims

    @contextmanager
    def open(self, **openparams):
        with self.signal.open(**openparams) as dset:
            yield dset

    @property
    def values(self):
        ret = np.empty(self.shape, dtype=self.dtype)
        with self.open(mode="r") as data:
            data.read_direct(ret)
        return ret


class NXRegularGrid(RegularGrid):
    def __init__(self, nxgroup, stackdim=None):
        """
        Args:
            nxgroup(nxfs.Path): NXdata or NXprocess
        """
        groups, axes, stackdim_prev = nxresult.regulargriddata(nxgroup)
        self.signals = [signal for group in groups.values() for signal in group]
        axnew = axis.factory(self.signals, type="nominal")
        if stackdim is None:
            stackdim = stackdim_prev
        axes.insert(stackdim, axnew)
        self.nxentry = nxgroup.nxentry()
        super(NXRegularGrid, self).__init__(axes, None, stackdim=stackdim)
        if stackdim <= stackdim_prev:
            stackdim_prev += 1
        self._stack_dims = [stackdim, stackdim_prev]

    @property
    def stack_dims(self):
        return self._stack_dims

    @contextmanager
    def open(self, **openparams):
        with self.nxentry.open(**openparams) as group:
            yield group

    @contextmanager
    def open_signals(self, **openparams):
        with ExitStack() as stack:
            yield [
                stack.enter_context(signal.open(**openparams))
                for signal in self.signals
            ]

    @property
    def signal_names(self):
        return [sig.name for sig in self.signals]

    def __getitem__(self, index):
        with self.open(mode="r") as group:

            def generator(signal):
                return group[self.nxentry.relpath(signal.path)]

            return indexing.getitem(
                generator,
                self.signals,
                index,
                self.ndim,
                shapefull=self.shape,
                axis=self.stackdim,
            )

    def __setitem__(self, index, value):
        with self.open(mode="a") as group:

            def selector(signal):
                return group[self.nxentry.relpath(signal.path)]

            indexing.setitem(
                selector,
                self.signals,
                index,
                self.ndim,
                value,
                shapefull=self.shape,
                axis=self.stackdim,
                method="set",
            )

    def __iadd__(self, value):
        with self.open() as group:

            def selector(signal):
                return group[self.nxentry.relpath(signal.path)]

            indexing.setitem(
                selector,
                self.signals,
                (Ellipsis,),
                self.ndim,
                value,
                shapefull=self.shape,
                axis=self.stackdim,
                method="add",
            )

    def __isub__(self, value):
        with self.open() as group:

            def selector(signal):
                return group[self.nxentry.relpath(signal.path)]

            indexing.setitem(
                selector,
                self.signals,
                (Ellipsis,),
                self.ndim,
                value,
                shapefull=self.shape,
                axis=self.stackdim,
                method="sub",
            )

    def __imul__(self, value):
        with self.open() as group:

            def selector(signal):
                return group[self.nxentry.relpath(signal.path)]

            indexing.setitem(
                selector,
                self.signals,
                (Ellipsis,),
                self.ndim,
                value,
                shapefull=self.shape,
                axis=self.stackdim,
                method="mul",
            )

    def __idiv__(self, value):
        with self.open() as group:

            def selector(signal):
                return group[self.nxentry.relpath(signal.path)]

            indexing.setitem(
                selector,
                self.signals,
                (Ellipsis,),
                self.ndim,
                value,
                shapefull=self.shape,
                axis=self.stackdim,
                method="div",
            )

    @property
    def dtype(self):
        with self.signals[0].open(mode="r") as dset:
            return dset.dtype

    @property
    def values(self):
        ret = np.empty(self.shape, dtype=self.dtype)
        for i, signal in enumerate(self.signals):
            with signal.open(mode="r") as dset:
                ret[i] = dset[()]
                # No longer works in h5py 3.1:
                # dset.read_direct(ret, dest_sel=(i, Ellipsis))
        return ret
