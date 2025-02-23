import numpy as np
import os
import fabio


class edfmemmap:
    """Access edf data with memmaps (cannot handle certain things like compression)"""

    def __init__(self, filename, mode="r"):
        # f = fabio.edfimage.EdfImage(filename)
        f = fabio.open(filename)

        self.dtype = f.bytecode
        self.shape = f.shape
        self.ndim = len(self.shape)
        offset = f._frames[f.currentframe].start
        self.mdata = np.memmap(
            filename, dtype=self.dtype, offset=offset, shape=self.shape, order="C"
        )

        if f.swap_needed():
            self.mdata.byteswap(True)

    @property
    def data(self):
        return self.mdata

    def __getitem__(self, index):
        return self.data[index]


class edfimage:
    """Access edf data with fabio"""

    def __init__(self, filename, mode="r"):
        try:
            self.f = fabio.open(filename)
        except Exception as e:
            raise OSError("Fabio cannot open file " + repr(filename)) from e
        self.ndim = 2

    @property
    def dtype(self):
        return self.f.bytecode

    @property
    def shape(self):
        return self.f.shape

    @property
    def data(self):
        return self.f.data

    def __getitem__(self, index):
        return self.data[index]

    @property
    def header(self):
        return self.f.header


def saveedf(filename, data, header, overwrite=False):
    exists = os.path.exists(filename)
    if exists:
        if overwrite:
            os.remove(filename)
        else:
            raise IOError("File exists (overwrite=False): {}".format(filename))
    fabio.edfimage.EdfImage(data=data, header=header).write(filename)
