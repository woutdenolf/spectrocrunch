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

from . import edf

import logging
logger = logging.getLogger(__name__)



        
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

class xiadata(object):

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

    def __init__(self,xiaconfig=None):
        if xiaconfig is None:
            self._xiaconfig = xiadict()
        else:
            self._xiaconfig = xiaconfig

    def skipdetectors(self,lst):
        self._xiaconfig.skipdetectors(lst)

    def keepdetectors(self,lst):
        self._xiaconfig.keepdetectors(lst)

    def onlydata(self):
        self._xiaconfig.onlydata()

    def onlystats(self):
        self._xiaconfig.onlystats()

    def dataandstats(self):
        self._xiaconfig.dataandstats()

    def overwrite(self,b):
        self._xiaconfig.overwrite(b)

    def onlyicrocr(self,b):
        self._xiaconfig.onlyicrocr(b)

    def dtcor(self,b):
        self._xiaconfig.dtcor(b)

    @property
    def nstats(self):
        if self._xiaconfig["indexing"]["onlyicrocr"]:
            return 2
        else:
            return self.NSTATS

    @property
    def data(self):
        self.onlydata()
        return self[:]

    @property
    def stats(self):
        self.onlystats()
        return self[:]

    def __getitem__(self, index):
        # Only data and no dtcor:  index applies on data
        # Only data and dtcor:  index applies on data, icrocr on stats
        # Only stats: index applies on stats
        #             index nonchanging: full or icrocr
        # Data and stats: index applies on data
        #                 index applies on stats but index[-2] must be full range

        if self._xiaconfig["indexing"]["data"]:
            # ...,nchan,ndet
            data = self._getdata(index)

            if not self._xiaconfig["indexing"]["stats"] and not self._xiaconfig["dtcor"]:
                return data

        if self._xiaconfig["indexing"]["stats"] or self._xiaconfig["dtcor"]:
            # ...,nstat,ndet
            if self._xiaconfig["indexing"]["data"]:
                index = indexing.replacefull(index,self.ndim,-2)
            else:
                if self._xiaconfig["indexing"]["onlyicrocr"]:
                    index = indexing.replace(index,self.ndim,-2,[self.STICR,self.STOCR])
                else:
                    index = indexing.replacefull(index,self.ndim,-2)

            stats = self._getstats(index)

            if not self._xiaconfig["indexing"]["data"]:
                return stats

        if self._xiaconfig["dtcor"]:
            index = (Ellipsis,self.STICR,slice(None))
            icr = stats[index]
            index = (Ellipsis,self.STOCR,slice(None))
            ocr = stats[index]

            if icr.size==1:
                data = data * self.CORTYPE(icr[0]) / ocr[0]
            else:
                s = list(stats.shape)
                s[-2] = 1
                icr = icr.reshape(s)
                ocr = ocr.reshape(s)
                data = data * np.asarray(icr,dtype=self.CORTYPE) / ocr

            if not self._xiaconfig["indexing"]["stats"]:
                return data

        return data,stats

    def _getdata(self,index=slice(None)):
        raise NotImplementedError("xiadata should not be instantiated, use one of the derived classes")

    def _getstats(self,index=slice(None)):
        raise NotImplementedError("xiadata should not be instantiated, use one of the derived classes")


class xialine(xiadata):

    def __init__(self,**kwargs):
        super(xialine,self).__init__(**kwargs)

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
                t1 = edf.edfmemmap(files[0]).dtype    
                t2 = self.stype
                t3 = self.CORTYPE
                return type(t1(0)*t2(0)*t3(0))
            else:
                return edf.edfmemmap(files[0]).dtype

    @property
    def stype(self):
        f = self.statfilename()
        if f is None:
            return self.XIASTYPE
        else:
            return edf.edfmemmap(f).dtype

    @property
    def dshape(self):
        # nspec x nchan x ndet
        files = self.datafilenames_used()
        n = len(files)
        if n==0:
            return ()
        else:
            s = edf.edfmemmap(files[0]).shape
            return s+(n,)

    @property
    def sshape(self):
        # nspec x nstat x ndet
        f = self.statfilename()
        if f is None:
            return ()
        else:
            s = edf.edfmemmap(f).shape
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
        generator = lambda f: edf.edfmemmap(f).data
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
            stats = edf.edfmemmap(f).data
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
        raise NotImplementedError("xialine should not be instantiated, use one of the derived classes")

    def statfilename(self):
        raise NotImplementedError("xialine should not be instantiated, use one of the derived classes")


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

        super(xialine_number,self).__init__(**kwargs)

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

        super(xialine_files,self).__init__(**kwargs)

    def __str__(self):
        return "xialine{}: {}".format(str(self.dshape),self._fileformat)
        
    def datafilenames(self):
        """
        Args:
            None
        Returns:
            list(str)
        """
        return self._datafiles

    def statfilename(self):
        """
        Args:
            None
        Returns:
            str or None
        """
        return self._statfile


class xiaimage(xiadata):

    def __init(self,**kwargs):
        super(xiaimage,self).__init__(**kwargs)

    def __str__(self):
        return "xiaimage{}".format(str(self.dshape))

    def _needlines(self,lines):
        if len(lines)==0:
            raise RuntimeError("No lines were defined")

    @property
    def ndim(self):
        return 3

    @property
    def dtype(self):
        lines = self._getlines()
        if len(lines)==0:
            return self.XIADTYPE
        else:
            return lines[0].dtype

    @property
    def stype(self):
        lines = self._getlines()
        if len(lines)==0:
            return self.XIASTYPE
        else:
            return lines[0].stype

    @property
    def dshape(self):
        # nlines x nspec x nchan x ndet
        lines = self._getlines()
        nlines = len(lines)
        if nlines==0:
            return ()
        else:
            return (nlines,) + lines[0].dshape

    @property
    def sshape(self):
        # nlines x nspec x nstat x ndet
        lines = self._getlines()
        nlines = len(lines)
        if len(lines)==0:
            return ()
        else:
            return (nlines,) + lines[0].sshape

    def _getdata(self,index=slice(None)):
        """
        Args:
            index (slice or num):
        Returns:
            array: nlines x nspec x nchan x ndet
        """
        lines = self._getlines()
        if len(lines)==0:
            np.empty((0,0,0,0),dtype=self.XIADTYPE)

        pindex,index = indexing.extract(index,4,0)

        if indexing.nonchanging(pindex):
            lines = lines[pindex]

        return np.stack([line[index] for line in lines],axis=0)

    def _getstats(self,index=slice(None)):
        """
        Args:
            index (slice or num):
        Returns:
            array: nlines x nspec x nstats x ndet
        """
        lines = self._getlines()
        if len(lines)==0:
            np.empty((0,0,0,0),dtype=self.XIADTYPE)

        if index is not Ellipsis:
            lines = lines.__getitem__(index)
        return np.stack([line.stats for line in lines],axis=0)

    def save(self,data,xialabels,stats=None):
        """
        Args:
            data(array): dimensions = nlines x nspec x nchan x ndet
            xialabels(list): number of labels = ndet
            st(Optional(array)): dimensions = nlines x nspec x nstats x ndet
        Returns:
            None
        """

        nlines, nspec, nchan, ndet = data.shape
        lines = self._getlines()
        if len(lines)!=nlines:
            raise RuntimeError("The map being saved has {} lines while {} lines were expected".format(nlines,len(lines)))
            
        for i in range(nlines):
            if stats is not None:
                stat = stats[i,...]
            else:
                stat = None
            lines[i].save(data[i,...],xialabels,stats=stat)

    def _getlines(self):
        return self._lines

    def skipdetectors(self,dets):
        """
        Args:
            dets(list): detector to skip when reading
        Returns:
            list: formated labels
        """
        for l in self._lines:
            l.skipdetectors(dets)

    def keepdetectors(self,dets):
        """
        Args:
            dets(list): detector to keep when reading
        Returns:
            list: formated labels
        """
        for l in self._lines:
            l.keepdetectors(dets)

    def datafilenames(self):
        """
        Args:
            None
        Returns:
            list(str)
        """
        files = []
        for l in self._lines:
            files += l.datafilenames()
        return files

    def statfilenames(self):
        """
        Args:
            None
        Returns:
            list(str)
        """
        files = []
        for l in self._lines:
            files.append(l.statfilename())
        return 

class xiaimage_linenumbers(xiaimage):
    """XIA image determined by its map number and line numbers
    """

    def __init__(self,path,radix,mapnum,linenums,**kwargs):
        """
        Args:
            path(list): path
            radix(list): radix
            mapnum(num): map number
            linenums(list(num)): line numbers

        Returns:
            None
        """
        self._lines = [xialine_number(path,radix,mapnum,linenum,**kwargs) for linenum in linenums]

        super(xiaimage_linenumbers,self).__init__(**kwargs)


class xiaimage_files(xiaimage):
    """XIA image determined by a list of files
    """

    def __init__(self,files,**kwargs):
        """
        Args:
            files(list(str)): data and stat file names

        Returns:
            None
        """

        files = xiagrouplines(files)
        self._lines = [xialine_files(f,**kwargs) for linenum,f in files.items()]

        super(xiaimage_files,self).__init__(**kwargs)


class xiaimage_number(xiaimage):
    """XIA image determined by its map number
    """

    def __init__(self,path,radix,mapnum,**kwargs):
        """
        Args:
            path(list): path
            radix(list): radix
            mapnum(num): map number

        Returns:
            None
        """
        self.path = path
        self.radix = radix
        self.mapnum = mapnum
        self.linekwargs = kwargs
        self._fileformat = os.path.join(path,xiaformat_map(radix,mapnum))
        self._lines = []

        super(xiaimage_number,self).__init__(**kwargs)

    def search(self):
        files = glob(self._fileformat.format("??","[0-9][0-9][0-9][0-9]*"))
        if len(files)==0:
            return
        files = xiagrouplines(files)
        n = len(files.values()[0])
        self._lines = [xialine_files(f,**self.linekwargs) for linenum,f in files.item() if len(f) == n]

    def _getlines(self):
        if len(self._lines)==0:
            self.search()
        return self._lines

    def save(self,data,xialabels,stats=None):
        """
        Args:
            data(array): dimensions = nlines x nspec x nchan x ndet
            xialabels(list): number of labels = ndet
            st(Optional(array)): dimensions = nlines x nspec x nstats x ndet
        Returns:
            None
        """

        if len(self._lines)==0:
            nlines, nspec, nchan, ndet = data.shape
            self._lines = [xialine_number(self.path,self.radix,self.mapnum,linenum,**self.linekwargs) for linenum in range(nlines)]

        super(xiaimage_number,self).save(data,xialabels,stats=stats)

