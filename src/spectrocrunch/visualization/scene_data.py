from ..io import spec
from ..io import fs
from ..utils import instance
from ..utils import units
from ..math import linop
from ..instruments import configuration
from ..process.edfregulargrid import EDFRegularGrid
from ..process.h5regulargrid import NXRegularGrid

import numpy as np
import collections


class Base(object):
    def __init__(
        self, axis0name=None, axis1name=None, linop0=None, linop1=None, **kwargs
    ):
        self.instrument = configuration.getinstrument(**kwargs)
        if axis0name is None:
            axis0name = self.instrument.imageaxes[0].upper()
        if axis1name is None:
            axis1name = self.instrument.imageaxes[1].upper()
        self.axis0name = axis0name
        self.axis1name = axis1name
        if not linop0:
            linop0 = 1, 0
        elif instance.isscalar(linop0):
            linop0 = 1, linop0
        self.linop0 = linop.LinearOperator(*linop0)
        if not linop1:
            linop1 = 1, 0
        elif instance.isscalar(linop1):
            linop1 = 1, linop1
        linop1 = linop.LinearOperator(*linop1)
        self.linop1 = linop1

    def motorquantity(self, values, motorname):
        return units.Quantity(values, self.instrument.units[motorname])


class PointBase(Base):
    def setcoordinates(self):
        coord0 = []
        coord1 = []
        lbls = []
        for positions, labels in self:
            p0 = sum(
                self.motorquantity(x, motname)
                for motname, x in zip(self.instrument.imagemotors, positions)
                if self.instrument.imageaxes[0] in motname
            )
            p1 = sum(
                self.motorquantity(x, motname)
                for motname, x in zip(self.instrument.imagemotors, positions)
                if self.instrument.imageaxes[1] in motname
            )

            # Append positions and labels
            coord0 += instance.asarray(p0).tolist()
            coord1 += instance.asarray(p1).tolist()
            lbls += labels
        if not coord0:
            raise RuntimeError("No Base found")
        self._coordinates0 = units.asqarray(coord0)
        self._coordinates1 = units.asqarray(coord1)
        self.labels = lbls

    @property
    def coordinates0(self):
        return self.linop0(self._coordinates0)

    @property
    def coordinates1(self):
        return self.linop1(self._coordinates1)


class ImageBase(Base):
    def __init__(self, grid, items=None, **kwargs):
        self.instrument = configuration.getinstrument(**kwargs)
        self.grid = grid
        self.set_items(items)
        axis0name, axis1name = [ax.title_nounits for ax in self._grid_image_axes]
        kwargs["axis0name"] = kwargs.get("kwargs", axis0name)
        kwargs["axis1name"] = kwargs.get("kwargs", axis1name)
        kwargs["instrument"] = self.instrument
        super(ImageBase, self).__init__(**kwargs)

    def set_items(self, items):
        idxgrid = [slice(None)] * self.grid.ndim
        stackdim = self.grid.stackdim
        stackaxis = self.grid.axes[stackdim]
        suffix = stackaxis.unit_suffix
        stack_dims = self.grid.stack_dims
        item_indices = []
        item_labels = []
        if items is None:
            for i, label in enumerate(stackaxis):
                idxgrid[stackdim] = i
                item_indices.append(tuple(idxgrid))
                item_labels.append(str(label) + suffix)
        else:
            for item in items:
                if not instance.isarray(item):
                    item = (item,)
                item = item + (-1,) * max(len(stack_dims) - len(item), 0)
                for i, v in zip(stack_dims, item):
                    j = self.grid.axes[i].locate(v, detectindex=True)
                    idxgrid[i] = j
                    if i == stackdim:
                        label = self.grid.axes[i][j]
                        if isinstance(label, fs.Path):
                            label = label.name
                item_indices.append(tuple(idxgrid))
                item_labels.append(str(label) + suffix)
        self.item_indices = item_indices
        self.item_labels = item_labels

    @property
    def _grid_image_axes(self):
        axes = [self.grid.axes[i] for i in self.grid.other_dims]
        if self.transpose:
            return axes[::-1]
        else:
            return axes

    @property
    def transpose(self):
        ax0name = self.grid.axes[self.grid.other_dims[0]].name
        if ax0name in self.instrument.imageaxes:
            return ax0name != self.instrument.imageaxes[0]
        else:
            return False

    @property
    def axis0values(self):
        return self.linop0(self._grid_image_axes[0].values)

    @property
    def axis1values(self):
        return self.linop1(self._grid_image_axes[1].values)

    def displaydata(self, index=None):
        """
        Args:
            index(Optional(list)): e.g. [5,None,1]

        Returns:
            data(array): e.g. nrow x ncol x 2
            channels: e.g. [0,None,1]
            labels: e.g. ["group5","group1"]
        """
        if index is None:
            it = zip(self.item_labels, self.item_indices)
            nimages = nout = len(self.item_indices)
        else:
            it = (
                (
                    (None, None)
                    if i is None
                    else (self.item_labels[i], self.item_indices[i])
                )
                for i in instance.asarray(index)
            )
            nout = len(index)
            nimages = nout - index.count(None)

        shape = self.grid.shape
        shape = [shape[i] for i in self.grid.other_dims] + [nimages]
        data = np.zeros(shape, dtype=self.grid.dtype)
        labels = [""] * nimages
        channels = [None] * nout
        iout = 0
        for itemidx, (label, idxgrid) in enumerate(it):
            if label is None:
                continue
            labels[iout] = label
            data[..., iout] = self.grid[idxgrid]
            channels[itemidx] = iout
            iout += 1

        if self.transpose:
            data = np.swapaxes(data, 0, 1)
        return data, channels, labels

    def interpolate(self, p0, p1):
        ind = list(range(len(p0)))
        data = self.grid.interpolate(None, p0, p1, asgrid=True, degree=1)
        data = data[:, ind, ind]
        result = collections.OrderedDict()
        for label, values in zip(self.grid.axes[0].values, data):
            result[label] = values
        return result


class EDFStack(ImageBase):
    def __init__(self, filenames, items, **kwargs):
        """
        Args:
            filename(str|list(str)): list of edf file names
            items(list(str))
        """
        if instance.isstring(filenames):
            filenames = [filenames]
        grid = EDFRegularGrid(filenames, instrument=kwargs.get("instrument", None))
        super(EDFStack, self).__init__(grid, items=items, **kwargs)


class NexusStack(ImageBase):
    def __init__(self, nxgroup, items, **kwargs):
        """
        Args:
            nxgroup(nxfs.Path): NXdata or NXprocess
            items(list(str))
        """
        grid = NXRegularGrid(nxgroup)
        super(NexusStack, self).__init__(grid, items=items, **kwargs)


class XanesSpec(PointBase):
    def __init__(self, filenames, specnumbers, labels=None, **kwargs):
        """
        Args:
            filenames(str|list(str)): list of spec file names
            specnumbers(list|list(list)): empty list of numbers => all xanes spectra
            labels(Optional(list|list(list))): uses the spec numbers by default
        """
        super(XanesSpec, self).__init__(**kwargs)

        if instance.isstring(filenames):
            filenames = [filenames]
        self.filenames = filenames
        self.specnumbers = self.listoflists(specnumbers)
        self.speclabels = self.listoflists(labels)
        self.setcoordinates()

    @staticmethod
    def listoflists(lst):
        if lst is None:
            return [[]]
        if not instance.isarray(lst):
            lst = [lst]
        if lst:
            if not instance.isarray(lst[0]):
                lst = [lst]
        else:
            lst = [lst]
        return lst

    def __iter__(self):
        for filename, numbers, labels in zip(
            self.filenames, self.specnumbers, self.speclabels
        ):
            # Get motor positions for each number
            f = spec.spec(filename)
            if not numbers:
                numbers = f.extractxanesginfo(
                    keepsum=True, sumingroups=False, keepindividual=False
                )
                numbers = [k[0] for k in numbers if len(k) == 1]
            if not numbers:
                continue
            positions = zip(
                *[f.getmotorvalues(nr, self.instrument.imagemotors) for nr in numbers]
            )
            if not labels:
                labels = numbers
            yield positions, labels
