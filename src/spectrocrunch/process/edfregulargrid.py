import os
from contextlib import contextmanager
import numpy as np

from . import axis
from .regulargrid import RegularGrid
from ..io import spec
from ..io.edf import edfimage
from ..instruments import configuration
from ..utils import units


class EDFRegularGrid(RegularGrid):
    def __init__(self, filenames, labels=None, **kwargs):
        instrument = configuration.getinstrument(**kwargs)
        parser = spec.edfheader_parser(
            units=instrument.units,
            compensationmotors=instrument.compensationmotors,
            axesnamemap=instrument.imagemotors,
            **instrument.edfheaderkeys,
        )
        header = edfimage(filenames[0]).header
        info = parser(header)
        axes = info["axes"]
        if len(axes) != 2:
            raise RuntimeError("Axes not specified in header")

        stackdim = 0
        if labels is None:
            labels = [".".join(os.path.basename(f).split(".")[:-1]) for f in filenames]
        axnew = axis.factory(labels, type="nominal")
        axes.insert(stackdim, axnew)

        self.filenames = filenames
        super(EDFRegularGrid, self).__init__(axes, None, stackdim=stackdim)

    @contextmanager
    def open(self, **openparams):
        yield np.stack(
            [edfimage(f, **openparams).data for f in self.filenames], axis=self.stackdim
        )

    @property
    def dtype(self):
        return edfimage(self.filenames[0]).dtype

    @property
    def values(self):
        with self.open(mode="r") as data:
            return data

    def __setitem__(self, index, value):
        raise NotImplementedError()


class XIARegularGrid(RegularGrid):
    def __init__(self, xiastack, stacklabel="energylabel", **kwargs):
        instrument = configuration.getinstrument(**kwargs)
        parser = spec.edfheader_parser(
            units=instrument.units,
            compensationmotors=instrument.compensationmotors,
            axesnamemap=instrument.imagemotors,
            **instrument.edfheaderkeys,
        )

        axes = None
        for i, xiaimage in enumerate(xiastack):
            header = xiaimage.header(source=instrument.metadata)
            info = parser(header)
            if axes:
                axes[0][i] = info[stacklabel]
            else:
                s = xiastack.dshape
                stackvalue = info[stacklabel].magnitude
                stackunit = info[stacklabel].units
                detectors = [int(det) for det in xiastack.detectors_used]
                axes = [
                    np.full(s[0], stackvalue),
                    info["axes"][0],
                    info["axes"][1],
                    axis.arange(s[3]),
                    axis.factory(detectors),
                ]
        axes[0] = axis.factory(units.Quantity(axes[0], stackunit))

        super(XIARegularGrid, self).__init__(axes, xiastack, stackdim=0)

    def __setitem__(self, index, value):
        raise NotImplementedError()
