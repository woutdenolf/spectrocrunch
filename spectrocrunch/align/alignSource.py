# -*- coding: utf-8 -*-
#
#   Copyright (C) 2015 European Synchrotron Radiation Facility, Grenoble, France
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

import os
import h5py
import numpy as np
import fabio

from .types import dataType

class alignSource(object):
    """Interface to data stacks with images as a function of energy, rotation angle, ...
    """

    def __init__(self,source,sublist,stackdim=None):
        """HDF5:
            source: file name or hdf5 object
            sublist: dataset names
           Files:
            source: base directory
            sublist: file names
           numpy array or external hdf5:
            source: list of ndarray or dataset
            sublist: None
        """

        if stackdim is None:
            stackdim = 2
        elif stackdim < 0 or stackdim > 2:
            raise ValueError("Stack dimension must be 0, 1 or 2.")

        self.stackdim = stackdim

        if isinstance(source,str):
            tmp, ext = os.path.splitext(source)
            if ext == '.h5' or ext == '.hdf5' or ext == '.nxs':
                self.handle = h5py.File(source, "r")
                self.datasets = [self.handle[name] for name in sublist]
                s = self.datasets[0].shape
                if not all(set.shape == s for set in self.datasets):
                    raise ValueError("Datasets don't have the same size.")
                s = self.initshape(s)
                self.sourcetype = dataType.h5
                self.nimages = s[stackdim]
                self.imgsize = tuple(np.delete(s,stackdim))
                self.nsets = len(self.datasets)

            elif ext == '':
                n = len(sublist[0])
                if not all(len(files) == n for files in sublist):
                    raise ValueError("Datasets don't have the same size.")
                self.datasets = [[os.path.join(source,f) for f in files] for files in sublist]

                self.sourcetype = dataType.singlefile
                self.nimages = n
                f = fabio.open(self.datasets[0][0])
                self.imgsize = (f.dim2,f.dim1)
                self.nsets = len(self.datasets)

            else:
                raise ValueError("Source type is not implemented.")
        elif isinstance(source,h5py.File) or isinstance(source,h5py.Group):
            self.handle = source
            self.datasets = [self.handle[name] for name in sublist]
            s = self.datasets[0].shape
            if not all(set.shape == s for set in self.datasets):
                raise ValueError("Datasets don't have the same size.")
            s = self.initshape(s)

            self.sourcetype = dataType.h5
            self.nimages = s[stackdim]
            self.imgsize = tuple(np.delete(s,stackdim))
            self.nsets = len(self.datasets)
        elif isinstance(source,list):
            if isinstance(source[0],(np.ndarray,h5py.Dataset)):
                
                s = self.initshape(source[0].shape)
                self.datasets = source

                self.sourcetype = dataType.nparray if isinstance(source[0],np.ndarray) else dataType.h5ext
                self.nimages = s[stackdim]
                self.imgsize = tuple(np.delete(s,stackdim))
                self.nsets = len(self.datasets)

            else:
                raise ValueError("Source type is not implemented.")
        else:
            raise ValueError("Source type is not implemented.")

    def initshape(self,s):
        if len(s) > 3:
            raise ValueError("Datasets should have 3 dimensions.")
        if len(s) < 3:
            self.stackdim = 2
            return s+(1,)*(3-len(s))
        return s

    def readimg(self,datasetindex,imageindex):
        if self.sourcetype==dataType.h5 or self.sourcetype==dataType.h5ext or self.sourcetype==dataType.nparray:
            #data = np.take(self.datasets[datasetindex],imageindex,axis=self.stackdim) # reads the entire file?!
            if self.stackdim==0:
                data = self.datasets[datasetindex][imageindex,...]
            elif self.stackdim==1:
                data = self.datasets[datasetindex][:,imageindex,:]
            else:
                data = self.datasets[datasetindex][...,imageindex]
        elif self.sourcetype==dataType.singlefile:
            data = fabio.open(self.datasets[datasetindex][imageindex]).data
        else:
            raise ValueError("Source type is not implemented.")

        if len(data.shape)==1:
            data = data[...,np.newaxis]

        return data
       
    @property 
    def dtype(self):
        if self.sourcetype==dataType.h5 or self.sourcetype==dataType.h5ext or self.sourcetype==dataType.nparray:
            return self.datasets[0].dtype
        elif self.sourcetype==dataType.singlefile:
            return fabio.open(self.datasets[0][0]).bytecode
        else:
            raise ValueError("Source type is not implemented.")

    def readimgas(self,datasetindex,imageindex,dtype):
        return self.readimg(datasetindex,imageindex).astype(dtype)
