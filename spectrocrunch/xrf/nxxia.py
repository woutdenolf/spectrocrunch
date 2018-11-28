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
from collections import defaultdict
from ..io import xiaedf
from ..io import nxfs
from ..io import spec
from ..instruments.configuration import getinstrument

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
                #self._create_appdef()
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
    def mca(self):
        self._xiaobject.dtcor(False)
        return self._xiaobject.data
        
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

    def _parse_xiadata(self):
        nxinstrument = self.nxinstrument()
        self._xiaobject.dtcor(False)
        mca = self.mca
        for i,detname in enumerate(self._xiaobject.detectors_used):
            nxdetector = nxinstrument.nxdetector('mca'+detname)
            nxdetector['counts'].mkfile(data=mca[...,i])
            nxdetector.updated()

    def _create_appdef(self):
        nxapplication = self.nxapplication
        for path in self.nxinstrument.iter_is_nxclass('NXdetector'):
            detname = path.name
            nxapplication[detname].link(path)
        nxapplication.updated()
