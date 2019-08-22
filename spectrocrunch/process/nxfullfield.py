# -*- coding: utf-8 -*-

import numpy as np
from . import nxprocess
from ..fullfield import import_id21


class Task(nxprocess.Task):

    DEFAULT_STACKDIM = 0

    def _parameters_defaults(self):
        super(Task, self)._parameters_defaults()
        self.required_parameters |= {
            'stackdim',
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
        }
        self.parameters['stackdim'] = self.parameters.get(
            'stackdim', self.DEFAULT_STACKDIM)

    def _execute(self):
        parameters = self.parameters
        darklib = import_id21.darklibrary(parameters)
        data, flat1, flat2, keyindices, stackaxes = import_id21.dataflatlibrary(
            parameters)

        # Save stack axes values
        positioners = self.temp_nxresults.positioners()
        for ax in stackaxes:
            positioners.add_axis(ax['name'], ax['data'])
        axes = [ax['name'] for ax in stackaxes]

        # Create NXdata groups for transmission and flat-field stacks
        nxdata = self.temp_nxresults.nxdata('detector0')
        signaldata = 'sample'
        signalflat1 = None
        signalflat2 = None
        if parameters['keepflat']:
            signalflat1 = 'flat1'
            if flat2 is not None:
                signalflat2 = 'flat2'
        needflat = parameters['keepflat'] or parameters['normalize']

        # Save/normalize image per image
        keys = list(data.keys())
        dtype = eval(parameters["dtype"])
        dim = [0]*3
        stackdim = parameters["stackdim"]
        imgdim = [i for i in range(3) if i != stackdim]
        dsetslice = [slice(None)]*3

        with nxdata.open() as group:
            for i, keyindex in enumerate(keyindices):
                dsetslice[stackdim] = i
                key = keys[keyindex]

                # Get data (DU/sec)
                img = import_id21.getnormalizedimage(
                    data[key], darklib, parameters)
                if needflat:
                    imgflat1 = import_id21.getnormalizedimage(
                        flat1[key], darklib, parameters)
                    if flat2 is not None:
                        imgflat2 = import_id21.getnormalizedimage(
                            flat2[key], darklib, parameters)

                # Normalize
                if parameters['normalize']:
                    if flat2 is None:
                        flat = imgflat1
                    else:
                        flat = (imgflat1+imgflat2)/2.
                    img = -np.log(img/flat)

                # Allocate memory
                if i == 0:
                    dim[imgdim[0]] = img.shape[0]
                    dim[imgdim[1]] = img.shape[1]
                    dim[stackdim] = len(data)

                    nxdata.add_signal(name=signaldata, shape=dim,
                                      chunks=True, dtype=dtype)
                    nxdata.set_axes(*axes)
                    dsetsample = group[signaldata]
                    if signalflat1:
                        nxdata.add_signal(name=signalflat1,
                                          shape=dim, chunks=True, dtype=dtype)
                        dsetflat1 = group[signalflat1]
                        if signalflat2:
                            nxdata.add_signal(
                                name=signalflat2, shape=dim, chunks=True, dtype=dtype)
                            dsetflat2 = group[signalflat2]

                # Save slice
                index = tuple(dsetslice)
                dsetsample[index] = img
                if signalflat1 is not None:
                    dsetflat1[index] = imgflat1
                    if signalflat2 is not None:
                        if len(flat2[key]) == 0:
                            dsetflat2[index] = dsetflat1[dsetslice]
                        else:
                            dsetflat2[index] = imgflat2

            # Save additional processing info
            info = self.temp_nxresults.nxcollection('info')
            if parameters["normalize"]:
                info["sample units"].mkfile(data="dimensionless")
            else:
                info["sample units"].mkfile(data="DU/sec")
            info["flat units"].mkfile(data="DU/sec")
            info["dark current"].mkfile(data="subtracted")
