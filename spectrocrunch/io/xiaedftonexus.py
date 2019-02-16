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

from . import xiaedf
from . import nxfs
from . import spec
from ..utils import units
from ..instruments.configuration import getinstrument

logger = logging.getLogger(__name__)


class Converter(object):
    
    def __init__(self, h5file=None, nxentry=None, qxrfgeometry=None, include_counters=None,
                 exclude_counters=None, diodeI0=None, diodeIt=None, dtype=np.float32,**parameters):
        self.dtype = dtype
        if nxentry:
            if not isinstance(nxentry,nxfs.Path):
                if h5file is None:
                    raise ValueError('Specify hdf5 file name')
                nxentry = nxfs.Path('/',h5file=h5file).nxentry(name=nxentry)
        else:
            if h5file is None:
                raise ValueError('Specify hdf5 file name')
            nxentry = nxfs.Path('/',h5file=h5file).new_nxentry()
        self._nxentry = nxentry
        self._qxrfgeometry = qxrfgeometry
        self._parameters = parameters
        if include_counters:
            self.includefuncs = include_counters
        else:
            self.includefuncs = []
        if exclude_counters:
            self.excludefuncs = exclude_counters
        else:
            self.excludefuncs = []
        if diodeI0:
            self._diodeI0 = diodeI0
        else:
            self._diodeI0 = 'I0'
        if diodeIt:
            self._diodeIt = diodeIt
        else:
            self._diodeIt = 'It'
        self._reset()
    
    def _reset(self):
        self._counter_paths = defaultdict(lambda: defaultdict(lambda: None))
        self._counter_names = None
        self._counter_to_interpretation = None
        self._interpretation_to_counter = None
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
                self._add_plots()
            finally:
                self._reset()
        return self.nxentry
        
    @property
    def nxentry(self):
        return self._nxentry
        
    @property
    def measurement(self):
        # NXdata or NXcollection
        return self.nxentry.measurement(nxclass='NXdata')

    @property
    def nxinstrument(self):
        nxinstrument = self.nxentry.nxinstrument()
        name = nxinstrument['name']
        if not name.exists:
            instrument = self.instrument
            name = name.mkfile(data=instrument.longname)
            name.update_stats(short_name=instrument.shortname)
        return nxinstrument
        
    @property
    def positioners(self):
        return self.nxentry.positioners()
    
    @property
    def nxmonochromator(self):
        return self.nxentry.nxmonochromator()
    
    @property
    def application(self):
        return self.nxentry.application('xrf',definition='NXxrf')
    
    @property
    def instrument(self):
        return getinstrument(**self._parameters)

    @property
    def counter_names(self):
        if self._counter_names is None:
            self._counter_names = self._xiaobject.counternames()
        return self._counter_names

    def _counter_interpretation(self):
        if self._interpretation_to_counter is None:
            ret = {}
            ctrdict = self.instrument.counterdict
            ret['icr'] = ctrdict.get('xrficr',None)
            ret['ocr'] = ctrdict.get('xrfocr',None)
            ret['lt'] = ctrdict.get('xrflt',None)
            ret['rt'] = ctrdict.get('xrfrt',None)
            ret['dt'] = ctrdict.get('xrfdt',None)
            ret['i0'] = ctrdict.get(self._diodeI0+'_counts',None)
            ret['it'] = ctrdict.get(self._diodeIt+'_counts',None)
            self._interpretation_to_counter = ret
            self._counter_to_interpretation = {v:k for k,v in ret.items() if v}

    @property
    def counter_to_interpretation(self):
        self._counter_interpretation()
        return self._counter_to_interpretation

    @property
    def interpretation_to_counter(self):
        self._counter_interpretation()
        return self._interpretation_to_counter

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
                         'energy':instrument.edfheaderkeys.get('energylabel',None),
                         'speclabel':instrument.edfheaderkeys.get('speclabel',None),
                         'slowlabel':instrument.edfheaderkeys.get('slowlabel',None),
                         'fastlabel':instrument.edfheaderkeys.get('fastlabel',None),
                         'units':instrument.units,
                         'compensationmotors':instrument.compensationmotors,
                         'axesnamemap':instrument.imagemotors}
            parser = spec.edfheader_parser(**parseinfo)
            if instrument.metadata=='counters':
                source = self.counter_names[0]
            else:
                source = None
            h = self._xiaobject.header(source=source)
            self._header = parser(h)
        return self._header
        
    @property
    def preset_time(self):
        return self.normalize_quantity(self.header['time'],'s')

    @property
    def energy(self):
        return self.normalize_quantity(self.header['energy'],'keV')

    def normalize_quantity(self,x,u):
        return units.astype(x.to(u),self.dtype)

    def add_measurement(self,measurement,name,data=None,units=None):
        if measurement.is_nxclass('NXdata'):
            return measurement.add_signal(name=name,data=data,units=units)
        else:
            return measurement[name].mkfile(data=data,units=units)

    def _parse_counters(self):
        measurement = self.measurement
        counters = self._xiaobject.counters
        interpretation = self.counter_to_interpretation
        for i,name in enumerate(self.counter_names):
            ctrname,detname = xiaedf.xianameparser.xiaparselabel(name)
            if self._skip_counter(ctrname):
                continue
            if detname:
                ctrpathname = ctrname+'_'+detname
            else:
                ctrpathname = ctrname
                detname = 'counters'
            ctrpath = self.add_measurement(measurement,ctrpathname,data=counters[...,i,0])
            interp = interpretation.get(ctrname,None)
            if interp:
                self._counter_paths[detname][interp] = ctrpath
        
        application = self.application
        for k in ['i0','it']:
            path = self._counter_paths['counters'][k]
            if path is not None:
                application[k].link(path)
        
        self.nxmonochromator.energy = self.energy
        measurement.updated()

    def _skip_counter(self,ctrname):
        if self.excludefuncs:
            bexclude = any(skip(ctrname) for skip in self.excludefuncs)
        else:
            bexclude = False
        if self.includefuncs:
            binclude = any(keep(ctrname) for keep in self.includefuncs)
        else:
            binclude = True
        return not binclude or bexclude

    def _parse_positioners(self):
        positioners = self.positioners
        axes = self.header['axes']
        for ax in axes:
            values = self.normalize_quantity(ax.values,ax.units)
            positioners.add_axis(ax.name,values,title=ax.title)
        
    def _parse_xiastats(self):
        measurement = self.measurement
        application = self.application
        self.stats = None
        preset_time = self.preset_time
        units = self.instrument.units
        
        if self._qxrfgeometry is not None:
            energy = self.energy.to('keV').magnitude
            op,preset_time = self._qxrfgeometry.I0op(energy,expotime=preset_time)
            application['i0_to_flux_factor'].mkfile(data=self.normalize_quantity(op.m,'Hz'))
            application['i0_to_flux_offset'].mkfile(data=self.normalize_quantity(op.b,'Hz'))
            op,preset_time = self._qxrfgeometry.Itop(energy,expotime=preset_time)
            application['it_to_flux_factor'].mkfile(data=self.normalize_quantity(op.m,'Hz'))
            application['it_to_flux_offset'].mkfile(data=self.normalize_quantity(op.b,'Hz'))
            distance = self.normalize_quantity(self._qxrfgeometry.getxrfdistance(),'cm')
            activearea = self.normalize_quantity(self._qxrfgeometry.getxrfactivearea(),'cm^2')
        
        for i,detname in enumerate(self._xiaobject.detectors_used):
            mcaname = 'mca'+detname
            counter_paths = self._counter_paths[detname]

            lt = counter_paths['lt']
            rt = counter_paths['rt']
            if lt is None or rt is None:
                events = self.stats[...,self._xiaobject.STEVT,i].astype(self.dtype)

            if lt is None:
                with np.errstate(divide='ignore', invalid='ignore'):
                    lt = events/self.stats[...,self._xiaobject.STICR,i].astype(self.dtype)
                lt[~np.isfinite(lt)] = 0
                lt = self.add_measurement(measurement,'live_time_'+detname,data=lt,units='s')
            else:
                u = units.get(lt.name,None)
                if u is not None:
                    lt.update_stats(units=u)
            
            if rt is None:
                with np.errstate(divide='ignore', invalid='ignore'):
                    rt = events/self.stats[...,self._xiaobject.STOCR,i].astype(self.dtype)
                rt[~np.isfinite(rt)] = np.nan
                rt = self.add_measurement(measurement,'elapsed_time_'+detname,data=rt,units='s')
            else:
                u = units.get(rt.name,None)
                if u is not None:
                    rt.update_stats(units=u)
                    
            nxgroup = application.nxcollection(mcaname)
            nxgroup['live_time'].link(lt)
            nxgroup['elapsed_time'].link(rt)
            if self._qxrfgeometry is not None:
                nxgroup['distance'].mkfile(data=distance[i])
                nxgroup['active area'].mkfile(data=activearea[i])
            nxgroup['preset_time'].mkfile(data=preset_time)
            
        measurement.updated()
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
        
        # Memory info message
        nbytes_mem *= ndet
        nbytes_data *= ndet
        nchunks *= ndet
        gb = 1024**3
        if nbytes_data>gb:
            u = 'GB'
        else:
            u = 'MB'
        nbytes_data = units.Quantity(nbytes_data,'bytes').to(u)
        if nbytes_mem>gb:
            u = 'GB'
        else:
            u = 'MB'
        nbytes_mem = units.Quantity(nbytes_mem,'bytes').to(u)
        msg = 'Copy {:~} of MCA data in {} chunks (available memory = {:~})'.format(nbytes_data,nchunks,nbytes_mem)
        
        return ind,msg

    def _copy_mca(self,nxdetector):
        copyindex,copymsg = self._copy_indexing()
        mcasum = None
        if len(copyindex)>2:
            logger.warning(copymsg)
            shape = self._xiaobject.dshape[:-1]
            dtype = self._xiaobject.dtype
            openparams = {'shape': shape,'dtype': dtype,'chunks': True}
            with nxdetector['data'].open(**openparams) as dset:
                h5file = dset.file
                #mcasum = np.empty(self._xiaobject.dshape[:-2],dtype=dtype)
                for i in range(len(copyindex)-1):
                    islice = (slice(copyindex[i],copyindex[i+1]),Ellipsis,0)
                    oslice = (slice(copyindex[i],copyindex[i+1]),Ellipsis)
                    mca = self._xiaobject[islice]
                    dset.write_direct(mca,dest_sel=oslice)
                    h5file.flush()
                    #mcasum[oslice] = mca.sum(axis=-1)
        else:
            logger.info(copymsg)
            mca = self._xiaobject.data[...,0]
            nxdetector['data'].mkfile(data=mca,chunks=True)
            #mcasum = mca.sum(axis=-1)

        nxdetector['data'].update_stats(interpretation='spectrum')
        nxdetector.updated()
        return nxdetector['data'],mcasum
        
    def _parse_xiadata(self):
        nxinstrument = self.nxinstrument
        measurement = self.measurement
        application = self.application
        self._xiaobject.dtcor(False)
        self._xiaobject.onlydata()
        for detname in self._xiaobject.iter_detectors():
            mcaname = 'mca'+detname
            # Save MCA data for 1 detector
            nxdetector = nxinstrument.nxdetector(mcaname)
            path,mcasum = self._copy_mca(nxdetector)
            # Add link in application subentry
            nxgroup = application[mcaname].mkdir()
            nxgroup[path.name].link(path)
            # Add MCA sum signal
            if mcasum is not None:
                path = self.add_measurement(measurement,'mca_sum',data=mcasum)
                if mcasum.ndim==2:
                    path.update_stats(interpretation='image')
        nxinstrument.updated()

    def _add_plots(self):
        measurement = self.measurement
        if measurement.is_nxclass('NXdata'):
            nxdata = measurement
        else:
            nxdata = self.nxentry.nxdata(nxfs.DEFAULT_PLOT_NAME)
        axes = [ax.name for ax in self.header['axes']]
        nxdata.set_axes(*axes)
        if nxdata!=measurement:
            for path in self.measurement:
                nxdata.add_signal(path=path)
        nxdata.mark_default()
