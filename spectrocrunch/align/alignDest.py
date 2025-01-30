# -*- coding: utf-8 -*-

import os
import h5py
import numpy as np
import fabio
import glob
import sys

from .types import dataType
from ..utils.listtools import move


class alignDest(object):
    """Interface for storing Elastix alignment results (list of image stacks)."""

    def __init__(self, dest, names, extension, stackdim=None, overwrite=False):
        """HDF5:
         dest: file name or hdf5 object
         names: dataset names
         extension: dataset name extension
        Files:
         dest: base directory
         names: file names
         extension: file name extension
        numpy array or external hdf5:
         dest: list of ndarray or dataset
         names: None
         extension: None
        """

        if stackdim is None:
            stackdim = 2
        elif stackdim < 0 or stackdim > 2:
            raise ValueError("Stack dimension must be 0, 1 or 2.")

        self.stackdim = stackdim
        self.overwrite = overwrite

        if isinstance(dest, str):
            tmp, ext = os.path.splitext(dest)
            if ext == ".h5" or ext == ".hdf5" or ext == ".nxs":
                self.desttype = dataType.h5
                self.handle = h5py.File(dest, "a")
                self.dir = "/"
                self.names = [s.strip(r"\/") for s in names]
                self.ext = extension
                # Datasets: self.dir+self.names[i]+self.ext
                self.datasets = []
            elif ext == "":
                self.desttype = dataType.singlefile
                self.dir = dest
                self.names = [s.strip(r"\/") for s in names]
                self.ext = extension
                self.destdir = ""
                self.format = ""
            else:
                raise ValueError("Destination type is not implemented.")
        elif isinstance(dest, h5py.File) or isinstance(dest, h5py.Group):
            self.desttype = dataType.h5
            self.handle = dest
            self.dir = "/"
            self.names = [s.strip(r"\/") for s in names]
            self.ext = extension
            # Datasets: self.dir+self.names[i]+self.ext
            self.datasets = []
        elif isinstance(dest, list):
            if isinstance(dest[0], np.ndarray):
                self.desttype = dataType.nparray
                # Copy to transfer ownership
                # self.datasets = dest
                self.datasets = [d.copy() for d in dest]
                for i in range(len(dest)):
                    dest[i] = self.datasets[i]
            elif isinstance(dest[0], h5py.Dataset):
                self.desttype = dataType.h5ext
                self.datasets = dest
            else:
                raise ValueError("Destination type is not implemented.")
        else:
            raise ValueError("Destination type is not implemented.")

    def destshape(self, nimages, imgsize):
        return tuple(move([imgsize[0], imgsize[1], nimages], 2, self.stackdim))

    def h5_composesetname(self, path, name):
        if path == "/":
            return path + name + self.ext
        else:
            return path + "/" + name + self.ext

    def h5_setnames(self, path):
        return [self.h5_composesetname(path, name) for name in self.names]

    def h5_datasetsexist(self, path):
        return not all(
            self.h5_composesetname(path, name) not in self.handle for name in self.names
        )

    def h5_nonexistingsetnames(self):
        num = 0
        path = self.dir
        while self.h5_datasetsexist(path):
            num += 1
            if self.dir == "/":
                path = self.dir + "entry" + str(num)
            else:
                path = self.dir + str(num)
        return self.h5_setnames(path)

    def h5_datasetscheck(self, shape, dtype):
        return all(
            self.handle[self.h5_composesetname(self.dir, name)].shape == shape
            and self.handle[name].dtype == dtype
            for name in self.names
        )

    def h5_createdatasets(self, setnames, shape, dtype):
        self.datasets = [
            self.handle.create_dataset(name, shape, dtype=dtype, chunks=True)
            for name in setnames
        ]

    def h5_prepare(self, shape, dtype):
        if self.overwrite:
            if self.h5_datasetsexist(self.dir):
                if self.h5_datasetscheck(shape, dtype):
                    setnames = self.h5_setnames(self.dir)
                    self.datasets = [self.handle[name] for name in setnames]
                else:
                    setnames = self.h5_nonexistingsetnames()
                    self.h5_createdatasets(setnames, shape, dtype)
            else:
                setnames = self.h5_setnames(self.dir)
                self.h5_createdatasets(setnames, shape, dtype)
        else:
            setnames = self.h5_nonexistingsetnames()
            self.h5_createdatasets(setnames, shape, dtype)

    def singlefile_composefilename(self, datasetindex, imageindex):
        return (
            os.path.join(self.destdir, self.names[datasetindex])
            + self.format % imageindex
            + self.ext
        )

    def singlefile_datasetsexist(self, path):
        return all(
            len(glob.glob(os.path.join(path, s) + "*" + self.ext)) != 0
            for s in self.names
        )

    def singlefile_nonexistingdestination(self):
        num = 0
        path = self.dir
        while self.singlefile_datasetsexist(path):
            num += 1
            path = self.dir + "_" + str(num)
        return path

    def singlefile_prepare(self, nimages):
        if self.overwrite:
            self.destdir = self.dir
        else:
            self.destdir = self.singlefile_nonexistingdestination()

        if not os.path.exists(self.destdir):
            os.makedirs(self.destdir)

        self.format = "_%%0%dd" % int(np.floor(np.log10(nimages)) + 1)

    def nplike_datasetscheck(self, shape, dtype):
        return all(arr.shape == shape and arr.dtype == dtype for arr in self.datasets)

    def nparray_prepare(self, shape, dtype):
        if not self.nplike_datasetscheck(shape, dtype):
            if self.overwrite:
                for d in self.datasets:
                    if d.dtype != dtype:
                        raise ValueError("Destination should have type %s." % dtype)
                    np.ndarray.resize(d, shape, refcheck=False)
            else:
                raise ValueError("Destination should have shape (%s,%s,%s)." % shape)

    def h5ext_prepare(self, shape, dtype):
        if not self.nplike_datasetscheck(shape, dtype):
            if self.overwrite:
                for d in self.datasets:
                    if d.dtype != dtype:
                        raise ValueError("Destination should have type %s." % dtype)
                    d.resize(shape)
            else:
                raise ValueError(
                    "Destination should have shape (%d,%d,%d) and type %s."
                    % (shape + (dtype,))
                )

    def prepare(self, nimages, imgsize, dtype):
        if self.desttype == dataType.h5:
            shape = self.destshape(nimages, imgsize)
            self.h5_prepare(shape, dtype)
        elif self.desttype == dataType.singlefile:
            self.singlefile_prepare(nimages)
        elif self.desttype == dataType.nparray:
            shape = self.destshape(nimages, imgsize)
            self.nparray_prepare(shape, dtype)
        elif self.desttype == dataType.h5ext:
            shape = self.destshape(nimages, imgsize)
            self.h5ext_prepare(shape, dtype)
        else:
            raise ValueError("Destination type is not implemented.")

    def singlefile_getfabiohandle(self, img):
        filetype = {".tif": "tif", ".tiff": "tiff", ".edf": "edf"}
        if self.ext in filetype:
            klass_name = filetype[self.ext] + "image"
        else:
            klass_name = "edfimage"
            self.ext = ".edf"
        module = sys.modules.get("fabio." + klass_name, None)
        klass = getattr(module, klass_name)
        return klass(data=img)

    def writeimg(self, img, datasetindex, imageindex):
        if (
            self.desttype == dataType.h5
            or self.desttype == dataType.nparray
            or self.desttype == dataType.h5ext
        ):
            index = tuple(
                move([slice(None), slice(None), imageindex], 2, self.stackdim)
            )
            self.datasets[datasetindex][index] = img
        elif self.desttype == dataType.singlefile:
            h = self.singlefile_getfabiohandle(img)
            filename = self.singlefile_composefilename(datasetindex, imageindex)
            h.write(filename)
        else:
            raise ValueError("Destination type is not implemented.")
