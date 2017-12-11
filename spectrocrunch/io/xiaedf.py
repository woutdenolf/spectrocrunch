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
import fabio
import numpy as np
from copy import copy
import itertools
import collections
import operator
import numbers
import re

from ..common import indexing
from ..common import listtools
from ..common import instance
from ..common import cache

from . import edf

import logging
logger = logging.getLogger(__name__)

import contextlib

XiaName = collections.namedtuple('XiaName', ['radix', 'mapnum', 'linenum', 'label', 'baselabel', 'detector'])

class XiaNameParser():

    def __init__(self,counters=None):
        defaultcounters = ["PUZ_arr","PUZ_xmap","xmap","arr","zap"]
        
        number = "[0-9]{4,}"
        fnumber = "[0-9]{{4,}}"
        
        xiafmt = "^(?P<radix>.+)_(?P<label>xia..)_(?P<mapnum>{})_0000_(?P<linenum>{}).edf$"
        ctrfmt = "^(?P<radix>.+)_(?P<label>{{}}.+)_(?P<mapnum>{})_0000.edf$".format(fnumber)

        self.xianames = [re.compile(xiafmt.format(number,number))]
        self.xianames += [re.compile(ctrfmt.format(ctr)) for ctr in defaultcounters]
        if counters is not None:
            self.xianames += [re.compile(ctrfmt.format(ctr)) for ctr in counters]

    def _parse(self,filename):
        filename = os.path.basename(filename)
        for xianame in self.xianames:
            m = xianame.match(filename)
            if m:
                return m.groupdict()
        return {}

    @staticmethod
    def xiaparselabel(label):
        if label.startswith("xia"):
            return "xia",label[3:]
        elif any(label.startswith(ctr) for ctr in ["PUZ_xmap","xmap","zap"]):
            tmp = label.split("_")
            return '_'.join(tmp[:-1]),tmp[-1]
        else:
            return label,""

    def parse(self,filename):
        """
        Args:
            filename(str): path/[radix]_[label]_[num]_0000_[linenumber].edf  ->  spectra
                           path/[radix]_[label]_[num]_0000.edf               ->  counters
        Returns
            dict
        """
        
        m = self._parse(filename)
        if m:
            radix = m["radix"]
            label = m["label"]
            mapnum = int(m["mapnum"])
            if "linenum" in m:
                linenum = int(m["linenum"])
            else:
                linenum = -1
            baselabel,detector = self.xiaparselabel(label)
            return XiaName(radix=radix,mapnum=mapnum,linenum=linenum,label=label,baselabel=baselabel,detector=detector)
        else:
            raise RuntimeError("{} does not have the proper XIA format".format(filename))

xianameparser = XiaNameParser()

def xiafilename(radix,mapnum,linenum,label):
    if linenum==-1:
        return "{}_{}_{:04d}_0000.edf".format(radix,label,mapnum)
    else:
        return "{}_{}_{:04d}_0000_{:04d}.edf".format(radix,label,mapnum,linenum)

def xiaformat_line(radix,mapnum,linenum):
    return "{}_{}_{:04d}_0000_{:04d}.edf".format(radix,'{}',mapnum,linenum)

def xiaformat_line(radix,mapnum,linenum):
    return "{}_{}_{:04d}_0000_{:04d}.edf".format(radix,'{}',mapnum,linenum)

def xiaformat_map(radix,mapnum):
    return "{}_{}_{:04d}_0000_{}.edf".format(radix,'{}',mapnum,'{}')

def xiaformat_ctr(radix,mapnum):
    return "{}_{}_{:04d}_0000.edf".format(radix,'{}',mapnum)
    
def xiaformat_radix(radix):
    return "{}_{}_{}_0000_{}.edf".format(radix,'{}','{}','{}')
    
def xiaformat_radixctr(radix):
    return "{}_{}_{}_0000.edf".format(radix,'{}','{}')

def xiaformat():
    return "{}_{}_{}_0000_{}.edf"

def ctrformat():
    return "{}_{}_{}_0000.edf"
    
def xiasortkey(filename):
    """Get key for sorting xia files:

        >>> files.sort(key = xiasortkey)

       Sorting is based on:

        [path]/[radix(1)]_[label(4)]_[num(2)]_0000_[linenumber(3)].edf

    Args:
        filename(str): 
    Returns:
        tuple: sort key
    """

    o = xianameparser.parse(filename)
    return o.radix,o.mapnum,o.linenum,o.label

def xiasearch(path,radix=None,mapnum=None,linenum=None,label=None,sort=True,ctrlabel=None,ctrs=True,onlyctrs=False):
    if radix is None:
        radix = "*"
    if mapnum is None:
        mapnum = "[0-9][0-9][0-9][0-9]*"
    if instance.isnumber(mapnum):
        mapnum = "{:04d}".format(mapnum)
    if onlyctrs:
        ctrs = True
        
    if ctrs:
        if ctrlabel is None:
            ctrlabel = '*'
            
        # Cannot be done unless using regex (or xianameparser.parse afterwards like done here)
        mask = os.path.join(path,ctrformat().format(radix,ctrlabel,mapnum))
        ctrfiles = glob(mask)
        
        ctrfiles = [f for f in ctrfiles if xianameparser.parse(f).linenum==-1]

        if onlyctrs:
            ctrfiles.sort(key=xiasortkey)
            return ctrfiles
            
    if label is None:
        label = "xia*"
    if linenum is None:
        linenum = "[0-9][0-9][0-9][0-9]*"
    if instance.isnumber(linenum):
        linenum = "{:04d}".format(linenum)
        
    mask = os.path.join(path,xiaformat().format(radix,label,mapnum,linenum))
    files = glob(mask)
    
    if ctrs:
        files += ctrfiles
    
    if sort:
        files.sort(key=xiasortkey)

    return files

def xiadetectorselect_tostring(detectors):
    """Convert detector identifier to string:
        0 -> "xia00"
        "xia00" -> "xia00"
        "S0" -> "xiaS0"
        "st" -> "xiast"
        "xiast" -> "xiast"

    Args:
        detectors(list): 
    Returns:
        list(str)
    """

    lst = copy(detectors)
    for i,s in enumerate(detectors):
        if instance.isnumber(s):
            lst[i] = "xia{:02}".format(s)
        elif instance.isstring(s):
            if not s.startswith("xia"):
                lst[i] = "xia{}".format(s)
    return lst

def xiadetectorselect_tonumber(detectors):
    """Convert detector identifier to number:
        0 -> 0
        "xia00" -> 0
        "S0" -> skip
        "st" -> skip
        "xiast" -> skip

    Args:
        detectors(list): 
    Returns:
        list(str)
    """
    
    lst = []

    for s in detectors:
        if instance.isstring(s):
            if s.startswith("xia"):
                s = s[3:]
            if s.isdigit():
                lst.append(int(s))
        elif instance.isnumber(s):
            lst.append(s)
            
    return lst
            
def xiadetectorselect_files(filenames,skipdetectors,keepdetectors):
    """Select xia detectors

    Args:
        filenames(list(str)): filenames
        skip(list): xia labels to be skipped (numbers or strings)
        keep(list): xia labels to be kept (numbers or strings)
    Returns:
        list(str)
    """
    
    if len(filenames)==0:
        return filenames

    skip = xiadetectorselect_tostring(skipdetectors)
    keep = xiadetectorselect_tostring(keepdetectors)

    if not skip and not keep:
        return filenames
    
    if not skip:
        valid = lambda x: x in keep
    elif not keep:
        valid = lambda x: x not in skip
    else:
        valid = lambda x: x not in skip and x in keep
        
    filenames = [f for f in filenames if valid(xianameparser.parse(f).label)]
    
    return filenames

def xiadetectorselect_numbers(detectors,skipdetectors,keepdetectors):
    """Select xia detectors

    Args:
        detectors(list(int)): 
        skip(list): xia labels to be skipped (numbers or strings)
        keep(list): xia labels to be kept (numbers or strings)
    Returns:
        list(str)
    """

    if len(detectors)==0:
        return detectors
        
    skip = xiadetectorselect_tonumber(skipdetectors)
    keep = xiadetectorselect_tonumber(keepdetectors)

    if not skip and not keep:
        return detectors
    
    if not skip:
        valid = lambda x: x in keep
    elif not keep:
        valid = lambda x: x not in skip
    else:
        valid = lambda x: x not in skip and x in keep
        
    detectors = [i for i in detectors if valid(i)]

    return detectors
    
def xiagroupkey(filename):
    """Group sorted files like this:

        [path]/[radix(1)]_[label]_[num(2)]_0000_[linenumber(3)].edf

    Args:
        filename(str): 
    Returns:
        tuple: sort key
    """

    o = xianameparser.parse(filename)
    return o.radix,o.mapnum,o.linenum,filename
    
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

def xiagroupdetectorskey(filename):
    """Group sorted files like this:

        [path]/[radix(1)]_[label]_[num(2)]_0000_[linenumber(3)].edf

    Args:
        filename(str): 
    Returns:
        tuple: sort key
    """

    o = xianameparser.parse(filename)
    return o.detector,o.baselabel,filename
    
def xiagroupdetectors(files):
    """
    Args:
        files(list(str)): unsorted file names

    Returns:
        OrderedDict(OrderedDict(list(str))): ret[detector][baselabel] = ["...xmap_x1c_00...","...xmap_x1c_00...",...]
    """
    
    files.sort(key=xiasortkey)

    keys = [xiagroupdetectorskey(f) for f in files]

    ret = collections.OrderedDict()
    
    for detector, v0 in itertools.groupby(keys, operator.itemgetter(0)):
        ret[detector] = collections.OrderedDict()
        for baselabel, v1 in itertools.groupby(v0, operator.itemgetter(1)):
            ret[detector][baselabel] = [v[-1] for v in v1]

    return ret

def xiagroupmaps2(files):
    """
    Args:
        files(list(str)): unsorted file names

    Returns:
        OrderedDict(OrderedDict(OrderedDict(list(str)))): ret[radix_mapnum][linenum] = ["...xia00...","...xia01...",...]
    """

    files.sort(key=xiasortkey)

    keys = [xiagroupkey(f) for f in files]

    ret = collections.OrderedDict()
    
    for radix, v0 in itertools.groupby(keys, operator.itemgetter(0)):
        for mapnum, v1 in itertools.groupby(v0, operator.itemgetter(1)):
            k = "{}_{}".format(radix,mapnum)
            ret[k] = collections.OrderedDict()
            for linenum, v2 in itertools.groupby(v1, operator.itemgetter(2)):
                ret[k][linenum] = [v[-1] for v in v2]

    return ret

def xiagroupmapskey(filename):
    """Group sorted files like this:

        [path]/[radix]_xia[label]_[num(1)]_0000_[linenumber(2)].edf

    Args:
        filename(str): 
    Returns:
        tuple: sort key
    """

    o = xianameparser.parse(filename)
    return o.mapnum,o.linenum,filename

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

    o = xianameparser.parse(filename)
    return o.linenum,filename

def xiagrouplines(files):
    """
    Args:
        files(list(str)): unsorted file names with the same radix and mapnumber

    Returns:
        OrderedDict(list(str)): ret[linenum] = ["...xia00...","...xia01...",...]
                                ret[-1] = ["...ctr1...","...ctr2...",...]
    """

    files.sort(key=xiasortkey)

    keys = [xiagrouplineskey(f) for f in files]

    ret = collections.OrderedDict()
    
    for linenum, v0 in itertools.groupby(keys, operator.itemgetter(0)):
        ret[linenum] = [v[-1] for v in v0]

    return ret

def xiagroupnxiafiles(files):
    if isinstance(files,list):
        return len([f for f in files if xianameparser.parse(f).baselabel=="xia"])
    else:
        return sum([xiagroupnxiafiles(f) for k,f in files.items()])

class xiadict(dict):

    def __init__(self,*args,**kwargs):
        super(xiadict,self).__init__(*args,**kwargs)
        if "indexing" not in self:
            self["indexing"] = {"data":True,"stats":False,"onlyicrocr":True}
        else:
            if "data" not in self["indexing"]:
                self["indexing"]["data"] = True
            if "stats" not in self["indexing"]:
                self["indexing"]["stats"] = False
            if "onlyicrocr" not in self["indexing"]:
                self["indexing"]["onlyicrocr"] = True

        if "overwrite" not in self:
            self["overwrite"] = False

        if "detectors" not in self:
            self["detectors"] = {"skip":[],"keep":[],"add":False}
        else:
            if "skip" not in self["detectors"]:
                self["detectors"]["skip"] = []
            if "keep" not in self["detectors"]:
                self["detectors"]["keep"] = []
            if "add" not in self["detectors"]:
                self["detectors"]["add"] = False

        if "dtcor" not in self:
            self["dtcor"] = False
            
        if "norm" not in self:
            self["norm"] = {"ctr":None,"func":lambda x:x}
        else:
            if "ctr" not in "norm":
                self["norm"]["ctr"] = None
            if "func" not in "norm":
                self["norm"]["func"] = lambda x:x
                
        if "maxlevel" not in self:
            self["maxlevel"] = 0

        if "levelinfo" not in self:
            self["levelinfo"] = {"dt":0,"sum":0,"norm":1,"ctr":1}
            
        if "envactive" not in self:
            self["envactive"] = False

        if "counter_reldir" not in self:
            self["counter_reldir"] = "." # relative to the path of the xia files

    def skipdetectors(self,lst):
        self["detectors"]["skip"] = lst

    def keepdetectors(self,lst):
        self["detectors"]["keep"] = lst

    def onlydata(self):
        self["indexing"]["data"] = True
        self["indexing"]["stats"] = False
        self["indexing"]["counters"] = False
        
    def onlystats(self):
        self["indexing"]["data"] = False
        self["indexing"]["stats"] = True
        self["indexing"]["counters"] = False
        
    def onlycounters(self):
        self["indexing"]["data"] = False
        self["indexing"]["stats"] = False
        self["indexing"]["counters"] = True
                  
    def dataandstats(self):
        self["indexing"]["data"] = True
        self["indexing"]["stats"] = True
        self["indexing"]["counters"] = False
    
    def dataall(self):
        self["indexing"]["data"] = True
        self["indexing"]["stats"] = True
        self["indexing"]["counters"] = True
        
    def overwrite(self,b):
        self["overwrite"] = b

    def onlyicrocr(self,b):
        self["indexing"]["onlyicrocr"] = b

    def dtcor(self,b):
        self["dtcor"] = b
        
    def norm(self,ctr,func=None):
        self["norm"]["ctr"] = ctr
        if callable(func):
            self["norm"]["func"] = func
        else:
            self["norm"]["func"] = lambda x:x
        
    def detectorsum(self,b):
        self["detectors"]["add"] = b

    def envactive(self,b):
        self["envactive"] = b
    
    def keep_getitem(self):
        return self["indexing"]["data"],self["indexing"]["stats"]

    def restore_getitem(self,data,stats):
        self["indexing"]["data"] = data
        self["indexing"]["stats"] = stats

    def maxlevel(self,level):
        self["maxlevel"] = max(self["maxlevel"],level)

    def resetmaxlevel(self,maxlevel):
        self["maxlevel"] = 0

    def setlevelinfo(self):
        if self["envactive"]:
            return

        if self["indexing"]["stats"]:
            # otherwise, stats will be retrieved twice
            self["levelinfo"]["dt"] = self["maxlevel"]
        else:
            # saves memory
            self["levelinfo"]["dt"] = 0

        # Sum detector right after dt correction
        self["levelinfo"]["sum"] = self["levelinfo"]["dt"]

        # Normalization info is only available on level 1
        self["levelinfo"]["norm"] = 1

        # Counter info is only available on level 1
        self["levelinfo"]["ctr"] = 1
        
    def counter_reldir(self,reldir):
        self["counter_reldir"] = reldir
        
class cachedict(dict):

    def __init__(self,*args,**kwargs):
        super(cachedict,self).__init__(*args,**kwargs)
        if "imagehandles" not in self:
            self["imagehandles"] = cache.LimitedSizeDict(size_limit=20)

        if "levels" not in self:
            self["levels"] = dict()

    def flush(self):
        self["imagehandles"] = cache.LimitedSizeDict(size_limit=20)
        self["levels"] = dict()

    def cachelevel(self,level):
        if level not in self["levels"]:
            self["levels"][level] = dict()
        return self["levels"][level]

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

    def __init__(self,level,xiaconfig=None,cache=None):
        if xiaconfig is None:
            self._xiaconfig = xiadict()
        else:
            self._xiaconfig = xiaconfig

        if cache is None:
            self._cache = cachedict()
        else:
            self._cache = cache

        self._level = level
        self._xiaconfig.maxlevel(level)
        
        self._normfunc = None

    def _update_kwargs(self,kwargs):
        kwargs["xiaconfig"] = self._xiaconfig
        kwargs["cache"] = self._cache

    def _getedfimage(self,filename,cache=False):
        # Don't do any caching for the moment
        try:
            return edf.edfmemmap(filename)
        except:
            return edf.edfimage(filename)
        
        if filename in self._cache["imagehandles"]:
            img = self._cache["imagehandles"][filename]
        else:
            try:
                img = edf.edfmemmap(filename)
            except:
                img = edf.edfimage(filename)
            if cache:
                self._cache["imagehandles"][filename] = img
        return img

    def _levelcache(self):
        return self._cache.cachelevel(self._level)

    def skipdetectors(self,lst):
        self._xiaconfig.skipdetectors(lst)

    def keepdetectors(self,lst):
        self._xiaconfig.keepdetectors(lst)

    def onlydata(self):
        self._xiaconfig.onlydata()

    def onlystats(self):
        self._xiaconfig.onlystats()

    def onlycounters(self):
        self._xiaconfig.onlycounters()
        
    @property
    def _lowestlevel(self):
        return self._level == 0

    @property
    def _highestlevel(self):
        return self._level == self._xiaconfig["maxlevel"]

    def _levelinfo(self):
        if self._highestlevel:   
            self._xiaconfig.setlevelinfo() 

        return {"dt":self._xiaconfig["dtcor"] and self._level == self._xiaconfig["levelinfo"]["dt"],\
                "sum":self._xiaconfig["detectors"]["add"] and self._level == self._xiaconfig["levelinfo"]["sum"],\
                "norm":self._xiaconfig["norm"]["ctr"] is not None and self._level == self._xiaconfig["levelinfo"]["norm"],\
                "counters":self._level == self._xiaconfig["levelinfo"]["ctr"]}

    @contextlib.contextmanager
    def env_onlydata(self):
        p = self._xiaconfig.keep_getitem()
        self.onlydata()
        self._xiaconfig.envactive(True)
        yield
        self._xiaconfig.restore_getitem(*p)

    @contextlib.contextmanager
    def env_onlystats(self):
        p = self._xiaconfig.keep_getitem()
        self.onlystats()
        self._xiaconfig.envactive(True)
        yield
        self._xiaconfig.restore_getitem(*p)

    @contextlib.contextmanager
    def env_onlycounters(self):
        p = self._xiaconfig.keep_getitem()
        self.onlycounters()
        self._xiaconfig.envactive(True)
        yield
        self._xiaconfig.restore_getitem(*p)
        
    @contextlib.contextmanager
    def env_cache(self):
        if self._highestlevel:
            self._cache.flush()
        yield

    def _tracemsg(self,msg):
        print("{}{}".format(" "*self._level*2,msg))
        
    def _dbgmsg(self,msg):
        return
        #if self._highestlevel:
        #if self._level >= 1:
        #    self._tracemsg(msg)

    def dataandstats(self):
        self._xiaconfig.dataandstats()

    def dataall(self):
        self._xiaconfig.dataall()
        
    def overwrite(self,b):
        self._xiaconfig.overwrite(b)

    def onlyicrocr(self,b):
        self._xiaconfig.onlyicrocr(b)

    def dtcor(self,b):
        self._xiaconfig.dtcor(b)

    def norm(self,ctr,func=None):
        return self.globalnorm(ctr,func=func)

    def globalnorm(self,ctr,func=None):
        self._xiaconfig.norm(ctr,func=func)
        self._normfunc = None
        
    def localnorm(self,ctr,func=None):
        self._xiaconfig.norm(ctr)
        self._normfunc = func

    def counter_reldir(self,reldir):
        self._xiaconfig.counter_reldir(reldir)
        
    def _getnormdata(self):
        raise NotImplementedError("This object does not have normalization information")
    
    @property
    def _normop(self):
        if self._normfunc is None:
            return self._xiaconfig["norm"]["func"]
        else:
            return self._normfunc
    
    def _getindexednormalizer(self,index):
        norm = self._getnormdata()
        norm = self._index_norm(norm,index)
        norm = self._normop(norm.astype(self.CORTYPE))
        return norm
        
    def detectorsum(self,b):
        self._xiaconfig.detectorsum(b)

    @property
    def data(self):
        with self.env_onlydata():
            return self[:]

    @property
    def stats(self):
        with self.env_onlystats():
            return self[:]

    @property
    def counters(self):
        with self.env_onlycounters():
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
        
    def _getaxis(self,i,stats=False,counters=False):
        """Axes (...,nchan,ndet) can transpose after indexing. Get the new axis position.

        Args:
            index(index):

        Returns:
            index
        """
        levelcache = self._levelcache()
        if stats:
            axes = levelcache["index"]["axesstats"]
        elif counters:
            axes = levelcache["index"]["axescounters"]
        else:
            axes = levelcache["index"]["axesdata"]
        if i<0:
            i += self.ndim
        
        if i not in axes:
            j = listtools.where(axes,lambda x:isinstance(x, list))
            if len(j)==0:
                axes = []
            else:
                axes = axes[j[0]]
        if i in axes:
            return axes.index(i)
        else:
            return None
            
    def _index_norm(self,norm,index):
        """Convert data indexing to norm indexing (i.e. remove channel and detector indexing) and apply to norm

        Args:
            index(index):

        Returns:
            index
        """
        
        norm = indexing.expanddims(norm,self.ndim)

        fullaxes = [-2,-1]
        indexfull,postindexfull,singletonindex = indexing.replacefull_transform(index,fullaxes,self.ndim)
        # Note: cannot levelcache this because singleton values in index may change

        norm = norm[indexfull]
        norm = postindexfull(norm)
        selind = [0]*len(singletonindex.selaxes)
        norm = singletonindex(norm,[selind])[0]
        return norm
        
    def _index_stats(self,index):
        return indexing.replacefull(index,self.ndim,[-2])

    def _index_counters(self,index):
        return indexing.replacefull(index,self.ndim,[-2,-1])
        
    def _cacheaxesorder(self,index):
        """Cache index and its derived information

        Args:
            index(index):

        Returns:
            array
        """
        levelcache = self._levelcache()
        if "index" not in levelcache:
            # Note: singleton values in index may change but this does not affect 
            levelcache["index"] = dict()
            levelcache["index"]["axesdata"],_ = indexing.axesorder_afterindexing(index,self.ndim)
            levelcache["index"]["axesstats"],_ = indexing.axesorder_afterindexing(self._index_stats(index),self.ndim)
            levelcache["index"]["axescounters"],_ = indexing.axesorder_afterindexing(self._index_counters(index),self.ndim)
            
    def extract_stat(self,stats,index,statind):
        """Extract specific statistics after indexing

        Args:
            stats(array): after indexing
            index(index):
            statind(list): statistics to be extracted

        Returns:
            list(array)
        """
        fullaxes = [-2]
        indexfull,postindexfull,singletonindex = indexing.replacefull_transform(index,fullaxes,self.ndim)
        # Note: cannot levelcache this because singleton values in index may change
        
        # indexfull already applied

        tmp = postindexfull(stats)
        
        nexpected = len(singletonindex.selaxes)
        if nexpected==1:
            selind = [[i] for i in statind]
        else:
            selind = [[i,0] for i in statind]

        tmp = singletonindex(tmp,selind)
  
        return tmp
        
    def extract_icrocr(self,stats,index):
        """Extract icr and ocr after indexing

        Args:
            stats(array): after indexing
            index(index):

        Returns:
            list(array)
        """
        return self.extract_stat(stats,index,[self.indexicr,self.indexocr])
        
    def sumdata(self,data,axis):
        """Sum data (after indexing) along axis and preserve dimensions

        Args:
            data(array): after indexing
            axis(num):

        Returns:
            array
        """
        
        s = list(data.shape)
        i = self._getaxis(-1)
        data = data.sum(axis=i)
        s[i] = 1
        return data.reshape(s)

    def __getitem__(self, index):
        # index applies on data
        # index[-2] extended applies on stats
        
        with self.env_cache():
            levelinfo = self._levelinfo()

            #self._dbgmsg("")
            #self._dbgmsg(self)
            #self._dbgmsg(self._xiaconfig)
            #self._dbgmsg("levelinfo: {}".format(levelinfo))
            #self._dbgmsg("index: {}".format(index))
            
            # What needs to be done on this level?
            correct = levelinfo["dt"] or levelinfo["sum"] or levelinfo["norm"]
            
            returndata = self._xiaconfig["indexing"]["data"]
            returnstats = self._xiaconfig["indexing"]["stats"]
            returncounters = self._xiaconfig["indexing"]["counters"]

            onlydata = not returnstats and not returncounters
            onlystats = not returndata and not returncounters
            onlycounters = not returndata and not returnstats
            
            needdata = returndata
            needstats = returnstats or (levelinfo["dt"] and self._xiaconfig["dtcor"])
            needcounters = returncounters
            
            # Cache index and its derived information
            self._cacheaxesorder(index)
            
            ### Get data from lower levels
            if needdata:
                with self.env_onlydata():
                    data = self._getdata(index)

                    #self._dbgmsg("data:")
                    #self._dbgmsg(self.dshape)
                    #self._dbgmsg(index)
                    #self._dbgmsg(data.shape)

                if onlydata and not correct:
                    #self._dbgmsg("return data\n")
                    return data

            ### Get statistics from lower levels
            if needstats:
                with self.env_onlystats():
                    levelcache = self._levelcache()
                    indexstats = self._index_stats(index)
                    stats = self._getstats(indexstats)

                    #self._dbgmsg("stats:")
                    #self._dbgmsg(self.sshape)
                    #self._dbgmsg(indexstats)
                    #self._dbgmsg(stats.shape)

                if onlystats:
                    #self._dbgmsg("return stats\n")
                    return stats

            ### Get counters from lower levels
            if needcounters:
                with self.env_onlycounters():
                    levelcache = self._levelcache()
                    indexcounters = self._index_counters(index)
                    counters = self._getcounters(indexcounters)

                    #self._dbgmsg("counters:")
                    #self._dbgmsg(self.cshape)
                    #self._dbgmsg(indexcounters)
                    #self._dbgmsg(counters.shape)

                if onlycounters:
                    #self._dbgmsg("return counters\n")
                    return counters
                    
            ### Apply corrections

            if levelinfo["dt"]:
                #self._dbgmsg("dt...")
                
                # Extract icr/ocr
                icr,ocr = self.extract_icrocr(stats,index)
                binvalid = np.logical_or(icr<=0,ocr<=0)
                if binvalid.any():
                    icr[binvalid] = 1
                    ocr[binvalid] = 1

                cor = icr.astype(self.CORTYPE) / ocr.astype(self.CORTYPE)
                
                # Apply deadtime correction
                data = data*cor
                #self._dbgmsg(data.shape)

            if levelinfo["sum"]:
                #self._dbgmsg("sum...")

                # Sum along the last axis (detectors)
                data = self.sumdata(data,-1)
                #self._dbgmsg(data.shape)

            if levelinfo["norm"]:
                #self._dbgmsg("norm...")
                
                # Get normalizer
                cor = self._getindexednormalizer(index)

                # Apply normalization
                data = data / cor.astype(self.CORTYPE)
                #self._dbgmsg(data.shape)
            
            # Return tuple
            if onlydata:
                #self._dbgmsg("return corrected data\n")
                return data
            if onlystats:
                #self._dbgmsg("return stats\n")
                return stats
            if onlycounters:
                #self._dbgmsg("return counters\n")
                return counters
                
            ret = tuple()
            if returndata:
                #self._dbgmsg("return corrected data\n")
                ret = ret + (data,)
            if returnstats:
                #self._dbgmsg("return stats\n")
                ret = ret + (stats,)
            if returncounters:
                #self._dbgmsg("return counters\n")
                ret = ret + (counters,)
            #self._dbgmsg("return data,stats\n")
            return ret

    def _getdata(self,index=slice(None)):
        raise NotImplementedError("xiadata should not be instantiated, use one of the derived classes")

    def _getstats(self,index=slice(None)):
        raise NotImplementedError("xiadata should not be instantiated, use one of the derived classes")

    def search(self):
        pass

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
        
    def xiadetectorselect_files(self,lst):
        return xiadetectorselect_files(lst,self._xiaconfig["detectors"]["skip"],self._xiaconfig["detectors"]["keep"])

    def xiadetectorselect_numbers(self,lst):
        return xiadetectorselect_numbers(lst,self._xiaconfig["detectors"]["skip"],self._xiaconfig["detectors"]["keep"])
    
    def datafilenames(self):
        raise NotImplementedError("xiadata should not be instantiated, use one of the derived classes")

    def datafilenames_used(self):
        files = self.datafilenames()
        if len(files)==0:
            return files

        return self.xiadetectorselect_files(files)
        
class xiacompound(xiadata):
    """Implements XIA data compounding (image, image stacks, ...)
    """

    def __init__(self,level,**kwargs):
        xiadata.__init__(self,level,**kwargs)
        
        self._items = []

    def __iter__(self):
        return iter(self.getitems())

    def getitems(self):
        if len(self._items)==0:
            self.search()
        return self._items

    @property
    def dtype(self):
        items = self.getitems()
        if len(items)==0:
            return self.XIADTYPE
        else:
            return items[0].dtype

    @property
    def stype(self):
        items = self.getitems()
        if len(items)==0:
            return self.XIASTYPE
        else:
            return items[0].stype

    @property
    def ctype(self):
        items = self.getitems()
        if len(items)==0:
            return self.XIADTYPE
        else:
            return items[0].ctype
            
    @property
    def dshape(self):
        # ... x nspec x nchan x ndet
        items = self.getitems()
        nlines = len(items)
        if nlines==0:
            return ()
        else:
            s = items[0].dshape
            return (nlines,) + s

    @property
    def sshape(self):
        # ... x nspec x nstat x ndet
        items = self.getitems()
        nlines = len(items)
        if nlines==0:
            return ()
        else:
            return (nlines,) + items[0].sshape

    @property
    def cshape(self):
        # ... x nrow x ncol x nctr x 1
        items = self.getitems()
        nimages = len(items)
        if nimages==0:
            return ()
        else:
            return (nimages,) + items[0].cshape

    def _getdata(self,index=slice(None)):
        """
        Args:
            index (slice or num):
        Returns:
            array: ... x nspec x nchan x ndet
        """
        items = self.getitems()
        if len(items)==0:
            np.empty((0,)*self.ndim,dtype=self.XIADTYPE)

        # Generate sliced stack
        generator = lambda line: line
        data = indexing.slicedstack(generator,items,index,self.ndim,shapefull=self.dshape,axis=0)
        
        if data.size==0:
            data = self._reshapeemptydata(data)
            
        return data
        
    def _reshapeemptydata(self,data):
        # slicedstack does not call the sublevels when the compound dimension is reduced to zero
        if self._xiaconfig["detectors"]["add"] and self._level > self._xiaconfig["levelinfo"]["sum"]:
            s = list(data.shape)
            s[self._getaxis(-1)] = 1
            data = data.reshape(s)
            
        return data

    def _getcounters(self,index=slice(None)):
        """
        Args:
            index (slice or num):
        Returns:
            array: ... x nspec x nctrs x 1
        """
        items = self.getitems()
        if len(items)==0:
            np.empty((0,)*self.ndim,dtype=self.XIADTYPE)

        # Generate sliced stack
        generator = lambda image: image
        return indexing.slicedstack(generator,items,index,self.ndim,shapefull=self.cshape,axis=0)
        
    def _getstats(self,index=slice(None)):
        """
        Args:
            index (slice or num):
        Returns:
            array: ... x nspec x nstats x ndet
        """
        items = self.getitems()
        if len(items)==0:
            np.empty((0,)*self.ndim,dtype=self.XIADTYPE)

        # Generate sliced stack
        generator = lambda line: line
        return indexing.slicedstack(generator,items,index,self.ndim,shapefull=self.sshape,axis=0)
        
    def save(self,data,xialabels,stats=None,ctrs=None,ctrnames=None,ctrheaders=None):
        """
        Args:
            data(array): dimensions = ... x nspec x nchan x ndet
            xialabels(list): number of labels = ndet
            stats(Optional(array)): dimensions = ... x nspec x nstats x ndet
            ctrs(Optional(array)): dimensions = ... nrow x ncol x ncounters
            ctrnames(Optional(list)):
            ctrheaders(Optional(array)): dimensions = ... nrow x ncol
        Returns:
            None
        """
        items = self.getitems()
        nitems = data.shape[0]
        if len(items)!=nitems:
            raise RuntimeError("The xia compound being saved has {} items while {} items were expected".format(nlines,len(items)))

        bctrsave = False
        
        for i in range(nitems):
            if stats is None:
                statsi = None
            else:
                statsi = stats[i,...]

            if isinstance(items[i],xiacompound):
                if ctrs is None:
                    ctrsi = None
                else:
                    ctrsi = ctrs[i,...]

                if ctrheaders is None:
                    ctrheaderi = None
                else:
                    ctrheaderi = ctrheaders[i,...]

                items[i].save(data[i,...],xialabels,stats=statsi,ctrs=ctrsi,ctrnames=ctrnames,ctrheaders=ctrheaderi)
            else:
                # self: xiaimage, item[i]: xialine
                bctrsave = ctrs is not None
                items[i].save(data[i,...],xialabels,stats=statsi)
                
        if bctrsave:
            ctrheaders = np.asscalar(np.array(ctrheaders)) # because it will be a 0-d array
            for k,v in zip(ctrnames,np.moveaxis(ctrs, -1, 0)):
                img = fabio.edfimage.EdfImage(data=v,header=ctrheaders)
                self._write(img,self._ctrformat.format(k))
                    
    def datafilenames(self):
        """
        Args:
            None
        Returns:
            list(str)
        """
        items = self.getitems()
        files = []
        for l in items:
            files += l.datafilenames()
        return files
        
    def ctrfilenames(self,ctr=None):
        """
        Args:
            ctr(Optional(str or list)): counter base name
        Returns:
            list(str)
        """
        items = self.getitems()
        files = []
        for l in items:
            files += l.ctrfilenames(ctr=ctr)
        return files

    def statfilenames(self):
        """
        Args:
            None
        Returns:
            list(str)
        """
        items = self.getitems()
        files = []
        for l in items:
            files += l.statfilenames()
        return files

    def counterbasenames(self):
        items = self.getitems()
        if items:
            return items[0].counterbasenames()
        else:
            return []

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
                t1 = self._getedfimage(files[0],cache=True).dtype    
                t2 = self.stype
                t3 = self.CORTYPE
                return type(t1(0)*t2(0)*t3(0))
            else:
                return self._getedfimage(files[0],cache=True).dtype

    @property
    def stype(self):
        f = self.statfilename()
        if f is None:
            return self.XIASTYPE
        else:
            return self._getedfimage(f,cache=True).dtype

    @property
    def dshape(self):
        # nspec x nchan x ndet
        files = self.datafilenames_used()
        n = len(files)
        if n==0:
            return ()
        else:
            s = self._getedfimage(files[0],cache=True).shape
            return s+(n,)

    @property
    def sshape(self):
        # nspec x nstat x ndet
        f = self.statfilename()
        if f is None:
            return ()
        else:
            s = self._getedfimage(f,cache=True).shape
            ndet = s[1]//self.NSTATS
            ndet = len(self.xiadetectorselect_numbers(range(ndet)))
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
        generator = lambda f: indexing.expanddims(self._getedfimage(f).data,self.ndim-1)
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
            stats = indexing.expanddims(self._getedfimage(f).data,self.ndim-1)
            s = stats.shape
            nspec = s[0]
            ndet = s[1]//self.NSTATS
            stats = stats.reshape((nspec,ndet,self.NSTATS))
            stats = np.swapaxes(stats,1,2)

            # Select stats of some or all detectors
            ind = self.xiadetectorselect_numbers(range(ndet))
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

    def save(self,data,xialabels,stats=None):
        """
        Args:
            data(array): dimensions = nspec x nchan x ndet
            xialabels(list): number of labels = ndet
            stats(Optional(array)): dimensions = nspec x nstats x ndet
        Returns:
            None
        """
        ndet = data.shape[-1]
        header = {"xdata":1,"xcorr":0,"xnb":1}
        for i in range(ndet):
            header["xdet"] = i
            img = fabio.edfimage.EdfImage(data=data[...,i],header=header)
            self._write(img,self._fileformat.format(xialabels[i]),makedir=i==0)

        if stats is not None:
            header = {"xdata":3,"xcorr":0,"xnb":ndet,"xdet":' '.join(str(i) for i in range(ndet))}

            nstats = stats.shape[1] * stats.shape[2]
            tmp = stats.reshape(( stats.shape[0], nstats ))
            ind = np.arange(nstats).reshape(stats.shape[1:3]).T.reshape(nstats)
            tmp = tmp[...,ind]

            img = fabio.edfimage.EdfImage(data=tmp,header=header)
            self._write(img,self._fileformat.format("xiast"))

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

class xiaimage(xiacompound):
    """Compounds XIA lines
    """
    def __init__(self,**kwargs):
        xiacompound.__init__(self,1,**kwargs)
        self._ctrfiles = []
        self._ctrformat = None

    def __str__(self):
        return "xiaimage{}".format(str(self.dshape))

    @property
    def ndim(self):
        return 4

    @property
    def ctype(self):
        files = self.ctrfilenames()
        if files:
            return self._getedfimage(files[0],cache=True).dtype
        else:
            return self.XIADTYPE  

    @property
    def cshape(self):
        # nrow x ncol x nctrs x 1
        files = self.ctrfilenames()
        if files:
            s = self._getedfimage(files[0],cache=True).shape
            return s+(len(files),1)
        else:
            return ()
            
    def _getcounters(self,index=slice(None)):
        """
        Args:
            index(Optional(slice)):
        Returns:
            array: nrow x ncol x nctrs x 1 (before indexing)
        """

        # Select all counters
        files = self.ctrfilenames()
        if len(files)==0:
            return np.empty((0,0,0,0),dtype=self.XIADTYPE)

        # Generate sliced stack
        generator = lambda f: indexing.expanddims(self._getedfimage(f).data,self.ndim-1)
        counters = indexing.slicedstack(generator,files,index,self.ndim,shapefull=self.cshape,axis=2)

        return counters
 
    def ctrfilenames(self,ctr=None):
        """
        Args:
            ctr(Optional(str or list)): counter base names
        Returns:
            list(str)
        """
        if len(self._ctrfiles)==0:
            self.search()
            
        if ctr is None:
            return self._ctrfiles
            
        if instance.isstring(ctr):
            ctr = [ctr]
        
        return list(listtools.flatten([f for c in ctr for f in self._ctrfiles if xianameparser.parse(f).baselabel==c]))
        
    def normfilename(self):
        """
        Args:
            None
        Returns:
            list(str)
        """
        ctr = self._xiaconfig["norm"]["ctr"]
        if ctr is None:
            return None
        
        files = self.ctrfilenames(ctr=ctr)
        if len(files)==1:
            return files[0]
        else:
            return None
        
    def _getnormdata(self):
        f = self.normfilename()
        if f is None:
            m = np.ones(self.dshape[0:2])
            logger.warning("No \"{}\" counter for {}".format(self._xiaconfig["norm"]["ctr"],self))
        else:
            m = self._getedfimage(f)

        return m
    
    def ctrsearch(self):
        path = os.path.abspath(os.path.join(self.path,self._xiaconfig["counter_reldir"]))
        self._ctrfiles = xiasearch(path,radix=self.radix,mapnum=self.mapnum,onlyctrs=True)
    
    def counterbasenames(self):
        return [xianameparser.parse(f).baselabel for f in self._ctrfiles]
        
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
        files = glob(self._fileformat.format("*"))
        files.sort(key = xiasortkey)

        valid = lambda x: "xia" in x and x!= "xiast"
        self._datafiles = [f for f in files if valid(xianameparser.parse(f).label)]

        valid = lambda x: x == "xiast"
        f = [f for f in files if valid(xianameparser.parse(f).label)]
        if len(f)==1:
            self._statfile = f[0]
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
        o = xianameparser.parse(files[0])
        radix,mapnum,linenum = o.radix,o.mapnum,o.linenum
        path = os.path.dirname(files[0])
        self._fileformat = os.path.join(path,xiaformat_line(radix,mapnum,linenum))

        # data file names
        valid = lambda x: x!="xiast" and "xia" in x
        self._datafiles = [f for f in files if valid(xianameparser.parse(f).label)]

        # stat file name
        tmp = [f for f in files if xianameparser.parse(f).label=='xiast']
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
        self._update_kwargs(kwargs)
        
        self.path = path
        self.radix = radix
        self.mapnum = mapnum
        self._ctrformat = os.path.join(self.path,xiaformat_ctr(self.radix,self.mapnum))
        
        self._items = [xialine_number(path,radix,mapnum,linenum,**kwargs) for linenum in linenums]

    def search(self):
        self.ctrsearch()
        
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
        self._update_kwargs(kwargs)

        if isinstance(files,list):
            files = xiagrouplines(files)

        self._items = [xialine_files(f,**kwargs) for linenum,f in files.items() if linenum != -1]
        
        if -1 in files:
            self._ctrfiles = files[-1]
            filename = self._ctrfiles[0]
        else:
            filename = files.values()[0][0]
        o = xianameparser.parse(filename) 
        self.radix, self.mapnum = o.radix, o.mapnum
        self.path = os.path.dirname(filename)
        self._ctrformat = os.path.join(self.path,xiaformat_ctr(self.radix,self.mapnum)) 
        
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
        self._update_kwargs(kwargs)

        self.path = path
        self.radix = radix
        self.mapnum = mapnum
        self._fileformat = os.path.join(self.path,xiaformat_map(self.radix,self.mapnum)) 
        self._ctrformat = os.path.join(self.path,xiaformat_ctr(self.radix,self.mapnum)) 
        
        self.linekwargs = kwargs

    def __str__(self):
        return "xiaimage{}: {}".format(str(self.dshape),self._fileformat)

    def search(self):
        files = xiasearch(self.path,radix=self.radix,mapnum=self.mapnum,ctrs=False)
        if len(files)!=0:
            files = xiagrouplines(files)
            n = xiagroupnxiafiles(files.values()[0]) # number of files in the first map
            self._items = [xialine_files(f,**self.linekwargs) for linenum,f in files.items() if xiagroupnxiafiles(f) == n]
        
        self.ctrsearch()
        
    def save(self,data,xialabels,stats=None,ctrs=None,ctrnames=None,ctrheaders=None):
        """
        Args:
            data(array): dimensions = nlines x nspec x nchan x ndet
            xialabels(list): number of labels = ndet
            stats(Optional(array)): dimensions = nlines x nspec x nstats x ndet
            ctrs(Optional(array)): dimensions = nlines x nspec x ncounters
            ctrnames(Optional(list)): 
            ctrheaders(Optional(array)): dimensions = ... nlines x nspec
        Returns:
            None
        """

        if len(self._items)==0:
            nlines = data.shape[0]
            self._items = [xialine_number(self.path,self.radix,self.mapnum,linenum,**self.linekwargs) for linenum in range(nlines)]
            if ctrnames is not None:
                self._ctrfiles = [self._ctrformat.format(k) for k in ctrnames]
                self._ctrfiles.sort(key=xiasortkey)
            
        xiaimage.save(self,data,xialabels,stats=stats,ctrs=ctrs,ctrnames=ctrnames,ctrheaders=ctrheaders)

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
        self._update_kwargs(kwargs)

        if isinstance(files,list):
            files = xiagroup(files)
        self._items = [xiaimage_files(fmap,**kwargs) for radix,fradix in files.items() for mapnum,fmap in fradix.items()]

class xiastack_radix(xiastack):
    """XIA stack determined by its radices
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
        self._update_kwargs(kwargs)

        if not isinstance(path,list):
            path = [path]
        if not isinstance(radix,list):
            radix = [radix]

        self.path = path
        self.radix = radix
        self._fileformat = [os.path.join(p,xiaformat_radix(r)) for p,r in zip(self.path,self.radix)]

        self.imagekwargs = kwargs

    def __str__(self):
        return "xiastack{}: {}".format(str(self.dshape),self._fileformat)

    def search(self):
        files = [xiasearch(p,radix=r) for p,r in zip(self.path,self.radix)]
        files = list(listtools.flatten(files))
        if len(files)==0:
            return
            
        files = xiagroupmaps2(files)

        n = xiagroupnxiafiles(files.values()[0]) # number of files in the first map
        self._items = [xiaimage_files(fmap,**self.imagekwargs) for radix_mapnum,fmap in files.items() if xiagroupnxiafiles(fmap)==n]
        
    def save(self,data,xialabels,stats=None,ctrs=None,ctrnames=None,ctrheaders=None):
        """
        Args:
            data(array): dimensions = nmaps x nlines x nspec x nchan x ndet
            xialabels(list): number of labels = ndet
            stats(Optional(array)): dimensions = nmaps x nlines x nspec x nstats x ndet
            ctrs(Optional(array)): dimensions = nmaps x nlines x nspec x ncounters
            ctrnames(Optional(list)): 
            ctrheaders(Optional(array)): dimensions = ... nlines x nspec
        Returns:
            None
        """

        if len(self._items)==0:
            nmaps = data.shape[0]
            self._items = [xiaimage_number(self.path[-1],self.radix[-1],mapnum,**self.imagekwargs) for mapnum in range(nmaps)]

        xiastack.save(self,data,xialabels,stats=stats,ctrs=ctrs,ctrnames=ctrnames,ctrheaders=ctrheaders)

class xiastack_mapnumbers(xiastack):
    """XIA stack determined by its radices and mapnumbers
    """

    def __init__(self,path,radix,mapnumbers,**kwargs):
        """
        Args:
            path(str or list(str)): path
            radix(str or list(str)): radix
            mapnumbers(list(int) or list(list(int))): mapnumbers

        Returns:
            None
        """
        xiastack.__init__(self,**kwargs)
        self._update_kwargs(kwargs)

        if not isinstance(path,list):
            path = [path]
        if not isinstance(radix,list):
            radix = [radix]
        if not isinstance(mapnumbers[0],list):
            mapnumbers = [mapnumbers]*len(radix) # each radix the same map numbers

        self.path = path
        self.radix = radix
        self.mapnumbers = dict(zip(radix,mapnumbers))
        
        self._fileformat = [os.path.join(p,xiaformat_radix(r)) for p,r in zip(self.path,self.radix)]

        self.imagekwargs = kwargs

    def __str__(self):
        return "xiastack{}: {}".format(str(self.dshape),self._fileformat)

    def search(self):
        files = [xiasearch(p,radix=r) for p,r in zip(self.path,self.radix)]
        files = list(listtools.flatten(files))
        if len(files)==0:
            return
            
        files = xiagroup(files)

        n = xiagroupnxiafiles(files.values()[0].values()[0]) # number of files in the first map
        self._items = [xiaimage_files(fmap,**self.imagekwargs) for radix,fradix in files.items()\
                                                     for mapnum,fmap in fradix.items()\
                                                     if mapnum in self.mapnumbers[radix] and xiagroupnxiafiles(fmap)==n]

    def save(self,data,xialabels,stats=None,ctrs=None,ctrnames=None,ctrheaders=None):
        """
        Args:
            data(array): dimensions = nmaps x nlines x nspec x nchan x ndet
            xialabels(list): number of labels = ndet
            stats(Optional(array)): dimensions = nmaps x nlines x nspec x nstats x ndet
            ctrs(Optional(array)): dimensions = nmaps x nlines x nspec x ncounters
            ctrnames(Optional(list)): 
            ctrheaders(Optional(array)): dimensions = ... nlines x nspec
        Returns:
            None
        """

        if len(self._items)==0:
            nmaps = data.shape[0]
            self._items = [xiaimage_number(self.path[-1],self.radix[-1],mapnum,**self.imagekwargs) for mapnum in range(nmaps)]

        xiastack.save(self,data,xialabels,stats=stats,ctrs=ctrs,ctrnames=ctrnames,ctrheaders=ctrheaders)
        
