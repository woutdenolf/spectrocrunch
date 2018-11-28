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

from __future__ import division
import numpy as np
from psutil import virtual_memory
from collections import defaultdict
import logging

from ..io import xiaedf
from ..io import nxfs
from ..io import spec
from ..utils import units
from ..instruments.configuration import getinstrument

logger = logging.getLogger(__name__)


def xiaimage_number_to_nx(path,radix,number,**kwargs):
    xiaimage = xiaedf.xiaimage_number(path,radix,number)
    convertor = Xia2NX(**kwargs)
    return convertor(xiaimage)


class Xia2NX(object):
    
    def __init__(self,h5file=None,nxentry=None,
                 counterstats=('icr','ocr','dt','lt','rt'),**parameters):
        if nxentry:
            if not isinstance(nxentry,nxfs.Path):
                if h5file is None:
                    raise ValueError('Specify hdf5 file name')
                nxentry = nxfs.Path('/',h5file=h5file).nxentry(name=entryname)
        else:
            if h5file is None:
                raise ValueError('Specify hdf5 file name')
            nxentry = nxfs.Path('/',h5file=h5file).new_nxentry()
        self._nxentry = nxentry
        self._parameters = parameters
        self._counterstats = counterstats
        self._reset()
    
    def _reset(self):
        self._counterinfo = defaultdict(lambda: defaultdict(lambda: None))
        self._counter_names = None
        self._header = None
        self._stats = None
        
    def __call__(self,xiaobject,**openparams):
        self._xiaobject = xiaobject
        with self.nxentry.open(**openparams):
            try:
                self._parse_counters()
                self._parse_positioners()
                self._parse_xiastats()
                self._parse_xiadata()
                self._create_appdef()
            finally:
                self._reset()
        return self.nxentry
        
    @property
    def nxentry(self):
        return self._nxentry
        
    @property
    def nxmeasurement(self):
        return self.nxentry.nxmeasurement()

    @property
    def nxinstrument(self):
        return self.nxentry.nxinstrument()
    
    @property
    def nxpositioners(self):
        return self.nxentry.positioners()
        
    @property
    def nxapplication(self):
        return self.nxentry.nxxrf('xrf')
    
    @property
    def instrument(self):
        return getinstrument(self._parameters)

    @property
    def counter_names(self):
        if self._counter_names is None:
            self._counter_names = self._xiaobject.counternames()
        return self._counter_names

    @property
    def stats(self):
        if self._stats is None:
            self._xiaobject.onlyicrocr(False)
            self._stats = self._xiaobject.stats
        return self._stats
    
    @stats.setter
    def stats(self,value):
        self._stats = value
    
    @property
    def header(self):
        if self._header is None:
            instrument = self.instrument
            parseinfo = {'time':instrument.edfheaderkeys.get('timelabel',None),
                         'speclabel':instrument.edfheaderkeys.get('speclabel',None),
                         'slowlabel':instrument.edfheaderkeys.get('slowlabel',None),
                         'fastlabel':instrument.edfheaderkeys.get('fastlabel',None),
                         'units':instrument.units}
            parser = spec.edfheader_parser(**parseinfo)
            h = self._xiaobject.header(source=self.counter_names[0])
            h = parser(h)
            if np.isnan(h['time'].magnitude):
                h = self._xiaobject.header()
                time = h['time']
            self._header = h
        return self._header
        
    @property
    def preset_time(self):
        return self.header['time'].to('s')

    def _parse_counters(self):
        nxmeasurement = self.nxmeasurement
        counters = self._xiaobject.counters
        for i,name in enumerate(self.counter_names):
            ctrpath = nxmeasurement[name].mkfile(data=counters[...,i,0])
            ctrname,detname = xiaedf.xianameparser.xiaparselabel(name)
            ctrname = ctrname.lower()
            for k in self._counterstats:
                if k in ctrname:
                    self._counterinfo[detname][k] = ctrpath
        nxmeasurement.updated()

    def _parse_positioners(self):
        nxpositioners = self.nxpositioners
        for ax in self.header['axes']:
            nxpositioners.add_axis(ax.name,ax.values,title=ax.title)

    def _parse_xiastats(self):
        nxinstrument = self.nxinstrument
        self.stats = None
        for i,detname in enumerate(self._xiaobject.detectors_used):
            detinfo = self._counterinfo[detname]
            nxdetector = nxinstrument.nxdetector('mca'+detname)

            dt = detinfo['dt']
            if dt is None:
                mask = self.stats[...,self._xiaobject.STDT,i]==100
            else:
                mask = dt.read()==100
            bmask = mask.any()
            
            events = self.stats[...,self._xiaobject.STEVT,i].astype(np.float32)

            lt = detinfo['lt']
            if lt is None:
                with np.errstate(divide='ignore', invalid='ignore'):
                    lt = events/self.stats[...,self._xiaobject.STICR,i]
            else:
                lt = lt.read()
            
            rt = detinfo['rt']
            if rt is None:
                with np.errstate(divide='ignore', invalid='ignore'):
                    rt = events/self.stats[...,self._xiaobject.STOCR,i]
            else:
                rt = rt.read()

            nxdetector['live_time'].mkfile(data=lt,units='s')
            nxdetector['elapsed_time'].mkfile(data=rt,units='s')
            nxdetector['preset_time'].mkfile(data=self.preset_time)
            nxdetector.updated()
        self.stats = None

    def _copy_indexing(self):
        # Index MCA array to avoid disk swapping
        shape = list(self._xiaobject.dshape)
        ndet = shape[-1]
        shape[-1] = self._xiaobject.dtype(0).itemsize
        nbytes_data = np.prod(shape)
        nbytes_mem = virtual_memory().available
        nchunks = int(np.ceil(nbytes_data/nbytes_mem))
        if nchunks>1:
            nchunks = int(np.ceil(2.*nbytes_data/nbytes_mem))
        n = shape[0]
        inc = int(np.ceil(n/nchunks))
        ind = range(n)[::inc]
        if ind[-1]!=n:
            ind.append(n)
        nbytes_data = units.Quantity(ndet*nbytes_data,'bytes').to('GB')
        nbytes_mem = units.Quantity(ndet*nbytes_mem,'bytes').to('GB')
        nchunks *= ndet
        msg = 'Copy {:~} of MCA data in {} chunks ({:~} available memory)'.format(nbytes_data,nchunks,nbytes_mem)
        return ind,msg

    def _copy_mca(self,nxdetector):
        copyindex,copymsg = self._copy_indexing()
        if len(copyindex)>2:
            logger.warning(copymsg)
            openparams = {'shape': self._xiaobject.dshape[:-1],
                          'dtype': self._xiaobject.dtype,
                          'chunks': True}
            with nxdetector['counts'].open(**openparams) as dset:
                for i in range(len(copyindex)-1):
                    islice = (slice(copyindex[i],copyindex[i+1]),Ellipsis)
                    mca = self._xiaobject[islice]
                    dset.write_direct(mca, source_sel=(Ellipsis,0),
                                           dest_sel=islice)
                    dset.file.flush()
                mca = self._xiaobject.data
                np.testing.assert_array_equal(dset[()],mca[...,0])
        else:
            logger.info(copymsg)
            mca = self._xiaobject.data
            nxdetector['counts'].mkfile(data=mca[...,0],chunks=True)
                
    def _parse_xiadata(self):
        nxinstrument = self.nxinstrument()
        self._xiaobject.dtcor(False)
        self._xiaobject.onlydata()
        for detname in self._xiaobject.iter_detectors():
            nxdetector = nxinstrument.nxdetector('mca'+detname)
            self._copy_mca(nxdetector)
            nxdetector.updated()

    def _create_appdef(self):
        nxapplication = self.nxapplication
        for path in self.nxinstrument.iter_is_nxclass('NXdetector'):
            detname = path.name
            nxapplication[detname].link(path)
        nxapplication.updated()
