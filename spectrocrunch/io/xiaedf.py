# -*- coding: utf-8 -*-
#
#   Copyright (C) 2017 European Synchrotron Radiation Facility, Grenoble, France
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
from glob import glob
import numbers
import fabio
import numpy as np
from copy import copy
import itertools
import collections
import operator
import numbers

from ..common import indexing
from ..common import listtools

from .edf import edfmemmap as edfimage

import logging
logger = logging.getLogger(__name__)

import contextlib


      
def xiaparsefilename(filename):
    """
    Args:
        filename(str): [path]/[radix]_xia[label]_[num]_0000_[linenumber].edf
    Returns
        tuple: radix, mapnum, linenum, label
    """
    lst = os.path.splitext(os.path.basename(filename))[0].split('_')
    return '_'.join(lst[:-4]),int(lst[-3]),int(lst[-1]),lst[-4][3:]

def xiafilename(radix,mapnum,linenum,label):
    return "{}_xia{}_{:04d}_0000_{:04d}.edf".format(radix,label,mapnum,linenum)

def xiaformat_line(radix,mapnum,linenum):
    return "{}_xia{}_{:04d}_0000_{:04d}.edf".format(radix,'{}',mapnum,linenum)

def xiaformat_map(radix,mapnum):
    return "{}_xia{}_{:04d}_0000_{}.edf".format(radix,'{}',mapnum,'{}')

def xiaformat_radix(radix):
    return "{}_xia{}_{}_0000_{}.edf".format(radix,'{}','{}','{}')

def xiaformat():
    return "{}_xia{}_{}_0000_{}.edf"

def xiasortkey(filename):
    """Get key for sorting xia files:

        >>> files.sort(key = xiasortkey)

       Sorting is based on:

        [path]/[radix(1)]_xia[label(4)]_[num(2)]_0000_[linenumber(3)].edf

    Args:
        filename(str): 
    Returns:
        tuple: sort key
    """

    radix,mapnum,linenum,label = xiaparsefilename(filename)

    if label.isdigit():
        label = int(label)

    return radix,mapnum,linenum,label

def xiasearch(path,radix=None,mapnum=None,linenum=None,label=None,sort=True):
    if radix is None:
        radix = '*'
    if mapnum is None:
        mapnum = "[0-9][0-9][0-9][0-9]*"
    if linenum is None:
        linenum = "[0-9][0-9][0-9][0-9]*"
    if label is None:
        label = "??"

    mask = os.path.join(path,xiaformat().format(radix,label,mapnum,linenum))

    files = glob(mask)
    if sort:
        files.sort(key=xiasortkey)

    return files

def xiadetectorselect(lst,skipdetectors,keepdetectors):
    """Select xia detectors

    Args:
        lst(list): 
        skip(list): xia labels to be skipped (numbers or strings)
        keep(list): xia labels to be kept (numbers or strings)
    Returns:
        list(str)
    """

    if len(lst)==0 or (len(skipdetectors)==0 and len(keepdetectors)==0):
        return lst

    if isinstance(lst[0],str):
        # Format: "00", "01", "S0", "S1", ...

        skip = copy(skipdetectors)
        for i in range(len(skip)):
            if isinstance(skip[i],numbers.Integral):
                skip[i] = "{:02}".format(skip[i])

        keep = copy(keepdetectors)
        for i in range(len(keep)):
            if isinstance(keep[i],numbers.Integral):
                keep[i] = "{:02}".format(keep[i])

        #Assume str format [path]/[radix]_xia[label]_[num]_0000_[linenumber].edf
        labels = [xiaparsefilename(f)[3] for f in lst]

        if len(skip)==0:
            lst = [f for l,f in zip(labels,lst) if l in keep]
        else:
            lst = [f for l,f in zip(labels,lst) if l not in skip or l in keep]

    else:
        # Format: 0, 1, ...

        skip = []
        for i in range(len(skipdetectors)):
            if isinstance(skipdetectors[i],str):
                if skipdetectors[i].isdigit(): #"00", "01", ... (there are no statistics for sums)
                    skip.append(int(skipdetectors[i]))
            else:
                skip.append(skipdetectors[i]) # number

        keep = []
        for i in range(len(keepdetectors)):
            if isinstance(keepdetectors[i],str):
                if keepdetectors[i].isdigit(): #"00", "01", ... (there are no statistics for sums)
                    keep.append(int(keepdetectors[i]))
            else:
                keep.append(keepdetectors[i]) # number

        if len(skip)==0:
            lst = [i for i in lst if i in keep]
        else:
            lst = [i for i in lst if i not in skip or i in keep]

    return lst

def xiagroupkey(filename):
    """Group sorted files like this:

        [path]/[radix(1)]_xia[label]_[num(2)]_0000_[linenumber(3)].edf

    Args:
        filename(str): 
    Returns:
        tuple: sort key
    """

    radix,mapnum,linenum,label = xiaparsefilename(filename)
    return radix,mapnum,linenum,filename

def xiagroup(files):
    """
    Args:
        files(list(str)): unsorted file names

    Returns:
        OrderedDict(OrderedDict(OrderedDict(list(str)))): ret[radix][mapnum][linenum] = ["...xia00...","...xia01...",...]
    """

    files.sort(key=xiasortkey)

    keys = [xiagroupkey(f) for f in files]

    ret = collections.OrderedDict()
    
    for radix, v0 in itertools.groupby(keys, operator.itemgetter(0)):
        ret[radix] = collections.OrderedDict()
        for mapnum, v1 in itertools.groupby(v0, operator.itemgetter(1)):
            ret[radix][mapnum] = collections.OrderedDict()
            for linenum, v2 in itertools.groupby(v1, operator.itemgetter(2)):
                ret[radix][mapnum][linenum] = [v[-1] for v in v2]

    return ret

def xiagroupmapskey(filename):
    """Group sorted files like this:

        [path]/[radix]_xia[label]_[num(1)]_0000_[linenumber(2)].edf

    Args:
        filename(str): 
    Returns:
        tuple: sort key
    """

    radix,mapnum,linenum,label = xiaparsefilename(filename)
    return mapnum,linenum,filename

def xiagroupmaps(files):
    """
    Args:
        files(list(str)): unsorted file names with the same radix

    Returns:
        OrderedDict(OrderedDict(list(str))): ret[mapnum][linenum] = ["...xia00...","...xia01...",...]
    """

    files.sort(key=xiasortkey)

    keys = [xiagroupmapskey(f) for f in files]

    ret = collections.OrderedDict()
    
    for mapnum, v0 in itertools.groupby(keys, operator.itemgetter(0)):
        ret[mapnum] = collections.OrderedDict()
        for linenum, v1 in itertools.groupby(v0, operator.itemgetter(1)):
            ret[mapnum][linenum] = [v[-1] for v in v1]

    return ret

def xiagrouplineskey(filename):
    """Group sorted files like this:

        [path]/[radix]_xia[label]_[num]_0000_[linenumber(1)].edf

    Args:
        filename(str): 
    Returns:
        tuple: sort key
    """

    radix,mapnum,linenum,label = xiaparsefilename(filename)
    return linenum,filename

def xiagrouplines(files):
    """
    Args:
        files(list(str)): unsorted file names with the same radix and mapnumber

    Returns:
        OrderedDict(list(str)): ret[linenum] = ["...xia00...","...xia01...",...]
    """

    files.sort(key=xiasortkey)

    keys = [xiagrouplineskey(f) for f in files]

    ret = collections.OrderedDict()
    
    for linenum, v0 in itertools.groupby(keys, operator.itemgetter(0)):
        ret[linenum] = [v[-1] for v in v0]

    return ret

def xiagroupnfiles(files):
    if isinstance(files,list):
        return len(files)
    else:
        return sum([xiagroupnfiles(f) for k,f in files.items()])

class xiadict(dict):

    def __init__(self,*args,**kwargs):
        super(xiadict,self).__init__(*args,**kwargs)
        if "indexing" not in self:
            self["indexing"] = {"data":True,"stats":False,"onlyicrocr":True}
        else:
            if "data" not in self["indexing"]:
                self["indexing"]["data"] = True
            if "stats" not in self["indexing"]:
                self["indexing"]["stats"] = True
            if "onlyicrocr" not in self["indexing"]:
                self["indexing"]["onlyicrocr"] = True

        if "overwrite" not in self:
            self["overwrite"] = False

        if "detectors" not in self:
            self["detectors"] = {"skip":[],"keep":[]}
        else:
            if "skip" not in self["detectors"]:
                self["detectors"]["skip"] = []
            if "keep" not in self["detectors"]:
                self["detectors"]["keep"] = []

        if "dtcor" not in self:
            self["dtcor"] = False

        if "maxlevel" not in self:
            self["maxlevel"] = 0

        if "dtcorlevel" not in self:
            self["dtcorlevel"] = 0

    def skipdetectors(self,lst):
        self["detectors"]["skip"] = lst

    def keepdetectors(self,lst):
        self["detectors"]["skip"] = lst

    def onlydata(self):
        self["indexing"]["data"] = True
        self["indexing"]["stats"] = False

    def onlystats(self):
        self["indexing"]["data"] = False
        self["indexing"]["stats"] = True

    def dataandstats(self):
        self["indexing"]["data"] = True
        self["indexing"]["stats"] = True

    def overwrite(self,b):
        self["overwrite"] = b

    def onlyicrocr(self,b):
        self["indexing"]["onlyicrocr"] = b

    def dtcor(self,b):
        self["dtcor"] = b

    def keep_getitem(self):
        return self["indexing"]["data"],self["indexing"]["stats"]

    def restore_getitem(self,data,stats):
        self["indexing"]["data"] = data
        self["indexing"]["stats"] = stats

    def maxlevel(self,level):
        self["maxlevel"] = max(self["maxlevel"],level)

    def resetmaxlevel(self,maxlevel):
        self["maxlevel"] = 0

    def setdtcorlevel(self):
        if self["indexing"]["stats"]:
            # otherwise, stats will be retrieved twice
            self["dtcorlevel"] = self["maxlevel"]
        else:
            # saves memory
            self["dtcorlevel"] = 0

class xiadata(object):
    """Implements XIA data access
    """

    STDET = 0
    STEVT = 1
    STICR = 2
    STOCR = 3
    STLT = 4
    STDT = 5
    NSTATS = 6
    XIADTYPE = np.int32
    XIASTYPE = np.int32
    CORTYPE = np.float64

    def __init__(self,level,xiaconfig=None):
        if xiaconfig is None:
            self._xiaconfig = xiadict()
        else:
            self._xiaconfig = xiaconfig

        self._level = level
        self._xiaconfig.maxlevel(level)

        self._cache_getitem = {"reducedstats":False,"axesstats":[],"axesdata":[],"channels_squeezed":False}

    def skipdetectors(self,lst):
        self._xiaconfig.skipdetectors(lst)

    def keepdetectors(self,lst):
        self._xiaconfig.keepdetectors(lst)

    def onlydata(self):
        self._xiaconfig.onlydata()

    def onlystats(self):
        self._xiaconfig.onlystats()

    @property
    def _lowestlevel(self):
        return self._level == 0

    @property
    def _highestlevel(self):
        return self._level == self._xiaconfig["maxlevel"]

    def _dtcor(self):
        if self._highestlevel:
            self._xiaconfig.setdtcorlevel()

        if self._xiaconfig["dtcor"]:
            dtcor = self._level == self._xiaconfig["dtcorlevel"]
        else:
            dtcor = False

        return dtcor

    @contextlib.contextmanager
    def env_onlydata(self):
        if self._highestlevel:
            p = self._xiaconfig.keep_getitem()
            self.onlydata()
        yield
        if self._highestlevel:
            self._xiaconfig.restore_getitem(*p)

    @contextlib.contextmanager
    def env_onlystats(self):
        if self._highestlevel:
            p = self._xiaconfig.keep_getitem()
            self.onlystats()
        yield
        if self._highestlevel:
            self._xiaconfig.restore_getitem(*p)

    def _dbgmsg(self,msg):
        return
        if self._highestlevel:
            print "{}{}".format(" "*self._level*2,msg)

    def dataandstats(self):
        self._xiaconfig.dataandstats()

    def overwrite(self,b):
        self._xiaconfig.overwrite(b)

    def onlyicrocr(self,b):
        self._xiaconfig.onlyicrocr(b)

    def dtcor(self,b):
        self._xiaconfig.dtcor(b)

    @property
    def data(self):
        self.onlydata()
        return self[:]

    @property
    def stats(self):
        self.onlystats()
        return self[:]

    @property
    def indexicr(self):
        if self._xiaconfig["indexing"]["onlyicrocr"]:
            return 0
        else:
            return self.STICR

    @property
    def indexocr(self):
        if self._xiaconfig["indexing"]["onlyicrocr"]:
            return 1
        else:
            return self.STOCR

    @property
    def nstats(self):
        if self._xiaconfig["indexing"]["onlyicrocr"]:
            return 2
        else:
            return self.NSTATS

    def _index_stats(self,index):
        """Convert data indexing to stat indexing (i.e. remove channel indexing)

        Args:
            index(index):

        Returns:
            index
        """
        return indexing.replacefull(index,self.ndim,-2)

    def _transpose_stats(self,index,stats):
        """Transpose stats after indexing in order to match data after indexing

        Args:
            index(index):
            stats(array): after indexing

        Returns:
            array
        """
        indexstats = self._index_stats(index)

        # Make stats match data after indexing:
        # dshape = ...,nchan,ndet
        # sshape = ...,nstats,ndet
        dimchannels = self.ndim-2
        axesdata,_ = indexing.axesorder_afterindexing(index,self.ndim)
        axesstats,_ = indexing.axesorder_afterindexing(indexstats,self.ndim)

        self._dbgmsg("axesdata {}".format(axesdata))
        self._dbgmsg("axesstats {}".format(axesstats))

        # Advanced indexing dimensions
        idatalist = listtools.where(axesdata,lambda x:isinstance(x, list))
        istatslist = listtools.where(axesstats,lambda x:isinstance(x, list))
        ndatalist = 0
        nstatlist = 0
        if len(idatalist)==0:
            ndatalist = 0
        else:
            ndatalist = len(axesdata[idatalist[0]])
        if len(istatslist)==0:
            nstatslist = 0
        else:
            nstatslist = len(axesstats[istatslist[0]])

        # Check where the channel dimension went
        channels_squeezed = dimchannels not in listtools.flatten(axesdata)
        
        if channels_squeezed:
            reshape = False

            # Put squeezed channel dimension at the end
            axesdata.append(dimchannels)
        else:
            reshape = True

            # Extract channel dimension when advanced
            if ndatalist!=0 and ndatalist!=nstatslist:
                
                if axesdata[idatalist[0]]==[dimchannels]:
                    reshape = True
                    axesdata[idatalist[0]] = dimchannels
                else:
                    reshape = False
                    axesdata[idatalist[0]] = [a for a in axesdata[idatalist[0]] if a!= dimchannels]
                    axesdata.append(dimchannels)

        # Transpose stats to match data after indexing
        ind = [axesstats.index(a) for a in axesdata if a in axesstats]
        self._dbgmsg("transpose {} ?".format(ind))
        if ind!=range(len(ind)) and sorted(ind)==range(len(ind)):
            self._dbgmsg("transpose {}".format(ind))
            stats = np.transpose(stats,ind)
            axesstats = listtools.listadvanced_int(axesstats,ind)

        self._dbgmsg("axesdata {}".format(axesdata))
        self._dbgmsg("axesstats {}".format(axesstats))

        # statistics index in stats for extracting
        statindex = axesstats.index(dimchannels)
        self._cache_getitem["statindex"] = statindex
        self._dbgmsg("statindex {}".format(statindex))
        
        if reshape:
            shape = list(stats.shape)
            shape[statindex] = 1
            shape = tuple(shape)
            self._dbgmsg("reshape {}".format(shape))
            self._cache_getitem["statreshape"] = lambda x:x.reshape(shape)
        else:
            self._cache_getitem["statreshape"] = lambda x:x

        return stats

    def extract_stat(self,stats,statind):
        """Extract specific statistics after indexing

        Args:
            stats(array): after indexing
            statind(list): specific statistics

        Returns:
            list(array)
        """

        i = self._cache_getitem["statindex"]
        reshape = self._cache_getitem["statreshape"]

        ndimstats = stats.ndim
        index = [slice(None)]*ndimstats

        ret = []
        for ind in statind:
            index[i] = ind
            stat = stats[index]
            stat = reshape(stat)
            ret.append(stat)

        return ret

    def extract_icrocr(self,stats):
        """Extract icr and ocr after indexing

        Args:
            stats(array): after indexing

        Returns:
            list(array)
        """
        return self.extract_stat(stats,[self.indexicr,self.indexocr])

    def __getitem__(self, index):
        # index applies on data
        # index[-2] extended applies on stats

        dtcor = self._dtcor()

        self._dbgmsg("")
        self._dbgmsg(self)
        self._dbgmsg(self._xiaconfig)
        self._dbgmsg("dtcor: {}".format(dtcor))

        if self._xiaconfig["indexing"]["data"]:
            with self.env_onlydata():
                data = self._getdata(index)

                self._dbgmsg("data:")
                self._dbgmsg(self.dshape)
                self._dbgmsg(index)
                self._dbgmsg(data.shape)

            if not self._xiaconfig["indexing"]["stats"] and not dtcor:
                self._dbgmsg("return data\n")
                return data

        if self._xiaconfig["indexing"]["stats"] or dtcor:
            with self.env_onlystats():
                indexstats = self._index_stats(index)
                stats = self._getstats(indexstats)

                self._dbgmsg("stats:")
                self._dbgmsg(self.sshape)
                self._dbgmsg(indexstats)
                self._dbgmsg(stats.shape)

                # Make stats match data after indexing:
                stats = self._transpose_stats(index,stats)

            if not self._xiaconfig["indexing"]["data"]:
                self._dbgmsg("return stats\n")
                return stats

        if dtcor:
            # Extract and reshape icr/ocr
            icr,ocr = self.extract_icrocr(stats)

            # Apply correction
            s = data.shape
            self._dbgmsg("data.shape {}".format(data.shape))
            self._dbgmsg("icr.shape {}".format(icr.shape))
            data = data * np.asarray(icr,dtype=self.CORTYPE) / np.asarray(ocr,dtype=self.CORTYPE)
            self._dbgmsg("data.shape {}".format(data.shape))
            data = data.reshape(s)

        if not self._xiaconfig["indexing"]["stats"]:
            self._dbgmsg("return dtcor data\n")
            return data

        self._dbgmsg("return data,stats\n")
        return data,stats

    def _getdata(self,index=slice(None)):
        raise NotImplementedError("xiadata should not be instantiated, use one of the derived classes")

    def _getstats(self,index=slice(None)):
        raise NotImplementedError("xiadata should not be instantiated, use one of the derived classes")


class xiacompound(xiadata):
    """Implements XIA data compounding (image, image stacks, ...)
    """

    def __init__(self,level,**kwargs):
        xiadata.__init__(self,level,**kwargs)

    def _getitems(self):
        if len(self._items)==0:
            self.search()
        return self._items

    def search(self):
        pass

    @property
    def dtype(self):
        items = self._getitems()
        if len(items)==0:
            return self.XIADTYPE
        else:
            return items[0].dtype

    @property
    def stype(self):
        items = self._getitems()
        if len(items)==0:
            return self.XIASTYPE
        else:
            return items[0].stype

    @property
    def dshape(self):
        # ... x nspec x nchan x ndet
        items = self._getitems()
        nlines = len(items)
        if nlines==0:
            return ()
        else:
            return (nlines,) + items[0].dshape

    @property
    def sshape(self):
        # ... x nspec x nstat x ndet
        items = self._getitems()
        nlines = len(items)
        if nlines==0:
            return ()
        else:
            return (nlines,) + items[0].sshape

    def _getdata(self,index=slice(None)):
        """
        Args:
            index (slice or num):
        Returns:
            array: ... x nspec x nchan x ndet
        """
        items = self._getitems()
        if len(items)==0:
            np.empty((0,)*self.ndim,dtype=self.XIADTYPE)

        # Generate sliced stack
        generator = lambda line: line
        return indexing.slicedstack(generator,items,index,self.ndim,shapefull=self.dshape,axis=0)

    def _getstats(self,index=slice(None)):
        """
        Args:
            index (slice or num):
        Returns:
            array: ... x nspec x nstats x ndet
        """
        items = self._getitems()
        if len(items)==0:
            np.empty((0,)*self.ndim,dtype=self.XIADTYPE)

        # Generate sliced stack
        generator = lambda line: line
        return indexing.slicedstack(generator,items,index,self.ndim,shapefull=self.sshape,axis=0)

    def save(self,data,xialabels,stats=None):
        """
        Args:
            data(array): dimensions = ... x nspec x nchan x ndet
            xialabels(list): number of labels = ndet
            st(Optional(array)): dimensions = ... x nspec x nstats x ndet
        Returns:
            None
        """
        items = self._getitems()
        nlines = data.shape[0]
        if len(items)!=nlines:
            raise RuntimeError("The xia compound being saved has {} items while {} items were expected".format(nlines,len(items)))
            
        for i in range(nlines):
            if stats is not None:
                stat = stats[i,...]
            else:
                stat = None
            items[i].save(data[i,...],xialabels,stats=stat)

    def datafilenames(self):
        """
        Args:
            None
        Returns:
            list(str)
        """
        items = self._getitems()
        files = []
        for l in items:
            files += l.datafilenames()
        return files

    def statfilenames(self):
        """
        Args:
            None
        Returns:
            list(str)
        """
        items = self._getitems()
        files = []
        for l in items:
            files += l.statfilenames()
        return files

class xialine(xiadata):

    def __init__(self,**kwargs):
        xiadata.__init__(self,0,**kwargs)

    def __str__(self):
        return "xialine{}".format(str(self.dshape))

    @property
    def ndim(self):
        return 3

    @property
    def dtype(self):
        files = self.datafilenames()
        if len(files)==0:
            if self._xiaconfig["dtcor"]:
                return self.CORTYPE
            else:
                return self.XIADTYPE
        else:
            if self._xiaconfig["dtcor"]:
                t1 = edfimage(files[0]).dtype    
                t2 = self.stype
                t3 = self.CORTYPE
                return type(t1(0)*t2(0)*t3(0))
            else:
                return edfimage(files[0]).dtype

    @property
    def stype(self):
        f = self.statfilename()
        if f is None:
            return self.XIASTYPE
        else:
            return edfimage(f).dtype

    @property
    def dshape(self):
        # nspec x nchan x ndet
        files = self.datafilenames_used()
        n = len(files)
        if n==0:
            return ()
        else:
            s = edfimage(files[0]).shape
            return s+(n,)

    @property
    def sshape(self):
        # nspec x nstat x ndet
        f = self.statfilename()
        if f is None:
            return ()
        else:
            s = edfimage(f).shape
            ndet = s[1]//self.NSTATS
            ndet = len(self.xiadetectorselect(range(ndet)))
            return (s[0],self.nstats,ndet)

    def _getdata(self,index=slice(None)):
        """
        Args:
            index(Optional(slice)):
        Returns:
            array: nspec x nchan x ndet (before indexing)
        """

        # Select some or all detectors
        files = self.datafilenames_used()
        if len(files)==0:
            return np.empty((0,0,0),dtype=self.XIADTYPE)

        # Generate sliced stack
        generator = lambda f: indexing.expanddims(edfimage(f).data,self.ndim-1)
        data = indexing.slicedstack(generator,files,index,self.ndim,shapefull=self.dshape,axis=-1)

        return data
    
    def _getstats(self,index=slice(None)):
        """
        Args:
            index(Optional(slice)):
        Returns:
            array: nspec x nstats x ndet (before indexing)
        """
        f = self.statfilename()
        if f is not None:
            # We have to read all stats
            stats = indexing.expanddims(edfimage(f).data,self.ndim-1)
            s = stats.shape
            nspec = s[0]
            ndet = s[1]//self.NSTATS
            stats = stats.reshape((nspec,self.NSTATS,ndet))

            # Select stats of some or all detectors
            ind = self.xiadetectorselect(range(ndet))
            if len(ind)==0:
                stats = np.empty((0,0,0),dtype=self.stype)
            elif len(ind)<ndet:
                stats = stats[...,ind]

            # Only ICR and OCR
            if self._xiaconfig["indexing"]["onlyicrocr"] and len(ind)!=0:
                stats = stats[:,[self.STICR,self.STOCR],:]

            # Apply index
            if not indexing.nonchanging(index):
                stats = stats[index]
        else:
            stats = np.empty((0,0,0),dtype=self.XIASTYPE)

        return stats

    def _write(self,img,filename,makedir=False):
        if not self._xiaconfig["overwrite"]:
            if os.path.isfile(filename):
                logger.warning("{} not saved (already exists)".format(filename))
                return

        if makedir:
            path = os.path.dirname(filename)
            if not os.path.exists(path):
                os.makedirs(path)

        img.write(filename)

    def save(self,data,xialabels,stats=None):
        """
        Args:
            data(array): dimensions = nspec x nchan x ndet
            xialabels(list): number of labels = ndet
            st(Optional(array)): dimensions = nspec x nstats x ndet
        Returns:
            None
        """
        for i in range(data.shape[-1]):
            img = fabio.edfimage.EdfImage(data=data[...,i])
            self._write(img,self._fileformat.format(xialabels[i]),makedir=i==0)

        if stats is not None:
            img = fabio.edfimage.EdfImage(data=stats.reshape(( stats.shape[0], stats.shape[1] * stats.shape[2] )))
            self._write(img,self._fileformat.format("st"))

    def xiadetectorselect(self,lst):
        return xiadetectorselect(lst,self._xiaconfig["detectors"]["skip"],self._xiaconfig["detectors"]["keep"])

    def datafilenames_used(self):
        files = self.datafilenames()
        if len(files)==0:
            return files

        return self.xiadetectorselect(files)

    def datafilenames(self):
        """
        Args:
            None
        Returns:
            list(str)
        """
        if len(self._datafiles)==0:
            self.search()
        return self._datafiles

    def statfilename(self):
        """
        Args:
            None
        Returns:
            str or None
        """
        if self._statfile is None:
            self.search()
        return self._statfile

    def statfilenames(self):
        return [self.statfilename()]

    def search(self):
        pass

class xiaimage(xiacompound):
    """Compounds XIA lines
    """
    def __init__(self,**kwargs):
        xiacompound.__init__(self,1,**kwargs)

    def __str__(self):
        return "xiaimage{}".format(str(self.dshape))

    @property
    def ndim(self):
        return 4

class xiastack(xiacompound):
    """Compounds XIA images
    """

    def __init__(self,**kwargs):
        xiacompound.__init__(self,2,**kwargs)

    def __str__(self):
        return "xiastack{}".format(str(self.dshape))

    @property
    def ndim(self):
        return 5

class xiamultistack(xiacompound):
    """Compounds XIA stacks
    """

    def __init__(self,**kwargs):
        xiacompound.__init__(self,3,**kwargs)

    def __str__(self):
        return "xiamultistack{}".format(str(self.dshape))

    @property
    def ndim(self):
        return 6

class xialine_number(xialine):
    """XIA line determined by its line number
    """
    
    def __init__(self,path,radix,mapnum,linenum,**kwargs):
        """
        Args:
            path(str): path
            radix(str): radix
            mapnum(num): map number
            linenum(num): line number
        Returns:
            None
        """

        self._fileformat = os.path.join(path,xiaformat_line(radix,mapnum,linenum))
        self._datafiles = []
        self._statfile = None

        xialine.__init__(self,**kwargs)

    def __str__(self):
        return "xialine{}: {}".format(str(self.dshape),self._fileformat)

    def search(self):
        """Find data and statistics files
        """
        files = glob(self._fileformat.format("[!s]*"))
        files.sort(key = xiasortkey)
        self._datafiles = files

        f = self._fileformat.format("st")
        if os.path.isfile(f):
            self._statfile = f
        else:
            self._statfile = None

class xialine_files(xialine):
    """XIA line determined by a list of files
    """

    def __init__(self,files,**kwargs):
        """
        Args:
            files(list(str)): data and stat files (order will be preserved)
        Returns:
            None
        """

        # filename format for saving
        radix,mapnum,linenum,_ = xiaparsefilename(files[0])
        path = os.path.dirname(files[0])
        self._fileformat = os.path.join(path,xiaformat_line(radix,mapnum,linenum))

        # data file names
        self._datafiles = [f for f in files if xiaparsefilename(f)[3]!='st']

        # stat file name
        tmp = [f for f in files if xiaparsefilename(f)[3]=='st']
        if len(tmp)!=0:
            self._statfile = tmp[0]
        else:
            self._statfile = None

        xialine.__init__(self,**kwargs)

    def __str__(self):
        return "xialine{}: {}".format(str(self.dshape),self._fileformat)

class xiaimage_linenumbers(xiaimage):
    """XIA image determined by its map number and line numbers
    """

    def __init__(self,path,radix,mapnum,linenums,**kwargs):
        """
        Args:
            path(str): path
            radix(str): radix
            mapnum(num): map number
            linenums(list(num)): line numbers

        Returns:
            None
        """

        xiaimage.__init__(self,**kwargs)
        kwargs["xiaconfig"] = self._xiaconfig

        self._items = [xialine_number(path,radix,mapnum,linenum,**kwargs) for linenum in linenums]

class xiaimage_files(xiaimage):
    """XIA image determined by a list of files
    """

    def __init__(self,files,**kwargs):
        """
        Args:
            files(list(str) or dict): data and stat file names

        Returns:
            None
        """

        xiaimage.__init__(self,**kwargs)
        kwargs["xiaconfig"] = self._xiaconfig

        if isinstance(files,list):
            files = xiagrouplines(files)

        self._items = [xialine_files(f,**kwargs) for linenum,f in files.items()]

class xiaimage_number(xiaimage):
    """XIA image determined by its map number
    """

    def __init__(self,path,radix,mapnum,**kwargs):
        """
        Args:
            path(str): path
            radix(str): radix
            mapnum(num): map number

        Returns:
            None
        """
        xiaimage.__init__(self,**kwargs)
        kwargs["xiaconfig"] = self._xiaconfig

        self.path = path
        self.radix = radix
        self.mapnum = mapnum
        self.linekwargs = kwargs
        self._fileformat = os.path.join(path,xiaformat_map(radix,mapnum))
        self._items = []

    def __str__(self):
        return "xiaimage{}: {}".format(str(self.dshape),self._fileformat)

    def search(self):
        files = glob(self._fileformat.format("??","[0-9][0-9][0-9][0-9]*"))
        if len(files)==0:
            return
        files = xiagrouplines(files)
        n = xiagroupnfiles(files.values()[0]) # number of files in the first map
        self._items = [xialine_files(f,**self.linekwargs) for linenum,f in files.items() if xiagroupnfiles(f) == n]

    def save(self,data,xialabels,stats=None):
        """
        Args:
            data(array): dimensions = nlines x nspec x nchan x ndet
            xialabels(list): number of labels = ndet
            st(Optional(array)): dimensions = nlines x nspec x nstats x ndet
        Returns:
            None
        """

        if len(self._items)==0:
            nlines = data.shape[0]
            self._items = [xialine_number(self.path,self.radix,self.mapnum,linenum,**self.linekwargs) for linenum in range(nlines)]

        xiaimage.save(self,data,xialabels,stats=stats)

class xiastack_files(xiastack):
    """XIA stack determined by a list of files
    """

    def __init__(self,files,**kwargs):
        """
        Args:
            files(list(str) or dict): data and stat file names

        Returns:
            None
        """

        xiastack.__init__(self,**kwargs)
        kwargs["xiaconfig"] = self._xiaconfig

        if isinstance(files,list):
            files = xiagroup(files)
        self._items = [xiaimage_files(fmap,**kwargs) for radix,fradix in files.items() for mapnum,fmap in fradix.items()]

class xiastack_radix(xiastack):
    """XIA stack determined by its radix
    """

    def __init__(self,path,radix,**kwargs):
        """
        Args:
            path(str or list(str)): path
            radix(str or list(str)): radix

        Returns:
            None
        """
        xiastack.__init__(self,**kwargs)
        kwargs["xiaconfig"] = self._xiaconfig

        if not isinstance(path,list):
            path = [path]
        if not isinstance(radix,list):
            radix = [radix]

        self.path = path
        self.radix = radix
        
        self._fileformat = [os.path.join(p,xiaformat_radix(r)) for p,r in zip(path,radix)]

        self.imagekwargs = kwargs
        self._items = []

    def __str__(self):
        return "xiastack{}: {}".format(str(self.dshape),self._fileformat)

    def search(self):
        files = [glob(f.format("??","[0-9][0-9][0-9][0-9]*","[0-9][0-9][0-9][0-9]*")) for f in self._fileformat]
        files = [f for f in listtools.flatten(files)]
        if len(files)==0:
            return
        files = xiagroupmaps(files)

        n = xiagroupnfiles(files.values()[0].values()[0]) # number of files in the first map
        self._items = [xiaimage_files(fmap,**self.imagekwargs) for radix,fradix in files.items() for mapnum,fmap in fradix.items() if xiagroupnfiles(fmap)==n]

    def save(self,data,xialabels,stats=None):
        """
        Args:
            data(array): dimensions = nmaps x nlines x nspec x nchan x ndet
            xialabels(list): number of labels = ndet
            st(Optional(array)): dimensions = nmaps x nlines x nspec x nstats x ndet
        Returns:
            None
        """

        if len(self._items)==0:
            nmaps = data.shape[0]
            self._items = [xiaimage_number(self.path[-1],self.radix[-1],mapnum,**self.imagekwargs) for mapnum in range(nmaps)]

        xiastack.save(self,data,xialabels,stats=stats)

