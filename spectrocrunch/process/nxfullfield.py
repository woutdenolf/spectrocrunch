# -*- coding: utf-8 -*-
#
#   Copyright (C) 2018 European Synchrotron Radiation Facility, Grenoble, France
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

import numpy as np
from . import nxtask
from ..io.nxfs import textarray
from ..fullfield import import_id21

class Task(nxtask.Task):
    
    DEFAULT_STACKDIM = 0
    
    def _parameters_defaults(self):
        super(Task,self)._parameters_defaults()
        self.allparams = [
                        # Input
                        "darklist",
                        "datalist",
                        "flatlist",
                        "flatbeforeafter",
                        # EDF header
                        "frametimelabel",
                        "frametimedefault",
                        "roilabel",
                        "nbframeslabel",
                        "stacklabel",
                        # Defaults
                        "dtype",
                        "darkcurrentzero",
                        "darkcurrentgain",
                        # Process
                        "normalize",
                        "keepflat",
                        "roi",
                        "rebin"
                        ]
        self._required_parameters(*self.allparams)
        self.parameters['stackdim'] = self.parameters.get('stackdim',self.DEFAULT_STACKDIM)
        
    def _parameters_filter(self):
        return super(Task,self)._parameters_filter()+self.allparams+['stackdim']
        
    def _execute(self):
        parameters = self.parameters
        darklib = import_id21.darklibrary(parameters)
        data,flat1,flat2,keyindices,stackaxes = import_id21.dataflatlibrary(parameters)

        # Save stack axes values
        positioners = self.nxresults.positioners()
        for ax in stackaxes:
            positioners.add_axis(ax['name'],ax['data'])
        axes = [ax['name'] for ax in stackaxes]
        
        # Create NXdata groups for transmission and flat-field stacks
        nxdata = self.nxresults.nxdata('detector0')
        signaldata = 'sample'
        signalflat1 = None
        signalflat2 = None
        if parameters['keepflat']:
            signalflat1 = 'flat1'
            if flat2 is not None:
                signalflat2 = 'flat2'
        needflat = parameters['keepflat'] or parameters['normalize']
        
        # Save/normalize image per image
        keys = data.keys()
        dtype = eval(parameters["dtype"])
        dim = [0]*3
        stackdim = parameters["stackdim"]
        imgdim = [i for i in range(3) if i!=stackdim]
        dsetslice = [slice(None)]*3
        
        with nxdata.open() as group:
            for i,keyindex in enumerate(keyindices):
                dsetslice[stackdim] = i
                key = keys[keyindex]
            
                # Get data (DU/sec)
                img = import_id21.getnormalizedimage(data[key],darklib,parameters)
                if needflat:
                    imgflat1 = import_id21.getnormalizedimage(flat1[key],darklib,parameters)
                    if flat2 is not None:
                        imgflat2 = import_id21.getnormalizedimage(flat2[key],darklib,parameters)
                
                # Normalize
                if parameters['normalize']:
                    if flat2 is None:
                        flat = imgflat1
                    else:
                        flat = (imgflat1+imgflat2)/2.
                    img = -np.log(img/flat)

                # Allocate memory
                if i==0:
                    dim[imgdim[0]] = img.shape[0]
                    dim[imgdim[1]] = img.shape[1]
                    dim[stackdim] = len(data)
                    
                    nxdata.add_signal(name=signaldata,shape=dim,chunks=True,dtype=dtype)
                    nxdata.set_axes(*axes)
                    dsetsample = group[signaldata]
                    if signalflat1:
                        nxdata.add_signal(name=signalflat1,shape=dim,chunks=True,dtype=dtype)
                        dsetflat1 = group[signalflat1]
                        if signalflat2:
                            nxdata.add_signal(name=signalflat2,shape=dim,chunks=True,dtype=dtype)
                            dsetflat2 = group[signalflat2]

                # Save slice
                index = tuple(dsetslice)
                dsetsample[index] = img
                if signalflat1 is not None:
                    dsetflat1[index] = imgflat1
                    if signalflat2 is not None:
                        if len(flat2[key])==0:
                            dsetflat2[index] = dsetflat1[dsetslice]
                        else:
                            dsetflat2[index] = imgflat2

            # Save additional processing info
            info = self.nxresults.nxcollection('info')
            if parameters["normalize"]:
                info["sample units"].mkfile(data=textarray("dimensionless"))
            else:
                info["sample units"].mkfile(data=textarray("DU/sec"))
            info["flat units"].mkfile(data=textarray("DU/sec"))
            info["dark current"].mkfile(data=textarray("subtracted"))
