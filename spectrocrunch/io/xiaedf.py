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

class xiadata(object):

    STDET = 0
    STEVT = 1
    STICR = 2
    STOCR = 3
    STLT = 4
    STDT = 5
    NSTATS = 6
    XIADTYPE = float
    XIASTYPE = int

    @property
    def data(self):
        return self._getdata()

    @property
    def stats(self):
        return self._getstats()

    @property
    def icrocr(self):
        stats = self._getstats()
        return stats[...,self.STICR,:],stats[...,self.STOCR,:]

    def _getdata(self,index=slice(None)):
        raise NotImplementedError("xiadata should not be instantiated, use one of the derived classes")

    def _getstats(self,index=slice(None)):
        raise NotImplementedError("xiadata should not be instantiated, use one of the derived classes")


class xialine(xiadata):

    def __init__(self,overwrite=False,skipdetectors=None,keepdetectors=None):
        """
        Args:
            path(list): path
            radix(list): radix
            num(numbers.Integral): map number
            linenum(numbers.Integral): line number
            linenum(Optional(overwrite)): line number
            skipdetectors(Optional(list)): detector numbers to skip when reading
            keepdetectors(Optional(list)): detector numbers to keep when reading
        Returns:
            None
        """

        self.overwrite = overwrite
        self.skipdetectors(skipdetectors)
        self.keepdetectors(keepdetectors)

    @property
    def dtype(self):
        files = self.datafilenames()
        if len(files)==0:
            return self.XIADTYPE
        else:
            return fabio.open(files[0]).bytecode

    @property
    def stype(self):
        f = self.statfilename()
        if f is None:
            return self.XIASTYPE
        else:
            return fabio.open(f).bytecode

    @property
    def dshape(self):
        # nspec x nchan x ndet
        files = self.datafilenames_used()
        n = len(files)
        if n==0:
            return ()
        else:
            s = fabio.open(files[0]).getDims()[::-1]
            s.append(n)
            return tuple(s)

    @property
    def sshape(self):
        # nspec x nstat x ndet
        f = self.statfilename()
        if f is None:
            return ()
        else:
            s = fabio.open(f).getDims()[::-1]
            ndet = s[1]//self.NSTATS
            ndet = len(self.xiadetectorselect(range(ndet)))
            return (s[0],self.NSTATS,ndet)

    def __getitem__(self, index):
        data = self._getdata(index)
        stats = self._getstats(indexing.expandindex(index,-2))
        return data,stats

    def skipdetectors(self,dets):
        """
        Args:
            dets(list): detector to skip when reading
        Returns:
            list: formated labels
        """
        if dets is None:
            self._skipdetectors = []
        else:
            self._skipdetectors = dets

    def keepdetectors(self,dets):
        """
        Args:
            dets(list): detector to keep when reading
        Returns:
            list: formated labels
        """
        if dets is None:
            self._keepdetectors = []
        else:
            self._keepdetectors = dets

    def _getdata(self,index=slice(None)):
        """
        Args:
            index(Optional(slice)):
        Returns:
            array: nspec x nchan x ndet
        """

        # python -mtimeit -s'import numpy as np' 'n=50' 'np.dstack([np.zeros((n,n))[0:n//2,0:n//2]]*n)'
        # python -mtimeit -s'import numpy as np' 'n=50' 'np.dstack([np.zeros((n,n))]*n)[0:n//2,0:n//2,:]'

        # Select some or all detectors
        files = self.datafilenames_used()
        if len(files)==0:
            return np.empty((0,0,0),dtype=self.XIADTYPE)

        if indexing.nonchangingindex(index):
            # Full range
            data = np.dstack([fabio.open(f).data for f in files])
        elif isinstance(index,tuple):
            # Separate indexing of detectors from indexing of the images
            findex = index[-1]
            if not indexing.nonchangingindex(findex):
                try:
                    files = files[findex]
                except:
                    files = [files[i] for i in findex]

                if len(files)==0:
                    return np.empty((0,0,0),dtype=self.XIADTYPE)
                elif isinstance(files,str):
                    files = [files]

            index = index[:-1]
            data = np.dstack([fabio.open(f).data[index] for f in files])
        else:
            data = np.dstack([fabio.open(f).data[index] for f in files])

        return data
    
    def _getstats(self,index=slice(None)):
        """
        Args:
            index(Optional(slice)):
        Returns:
            array: nspec x nstats x ndet
        """
        f = self.statfilename()
        if f is not None:
            # We have to read all stats
            stats = fabio.open(f).data
            s = stats.shape
            nspec = s[0]
            ndet = s[1]//self.NSTATS
            stats = stats.reshape((nspec,self.NSTATS,ndet))

            # Select stats of some or all detectors
            ind = self.xiadetectorselect(range(ndet))
            if len(ind)==0:
                stats = np.empty((0,0,0),dtype=self.dtype)
            elif len(ind)<ndet:
                stats = stats[...,ind]
            
            # Apply index
            if not indexing.nonchangingindex(index):
                stats = stats[index]

        else:
            stats = np.empty((0,0,0),dtype=self.XIASTYPE)

        return stats

    def _write(self,img,filename,makedir=False):
        if not self.overwrite:
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
        return xiadetectorselect(lst,self._skipdetectors,self._keepdetectors)

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

    def _needlines(self,lines):
        if len(lines)==0:
            raise RuntimeError("No lines were defined")

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

    def __getitem__(self, index):
        # index: tuple(slice) or slice or number

        if isinstance(index,tuple):
            pindex = index[0]
        else:
            pindex = index

        data = self._getdata(pindex)
        stats = self._getstats(pindex)

        if isinstance(index,tuple):
            data = data.__getitem__(index[1:])
            if stats is not None:
                stats = stats.__getitem__(index[1:])

        return data,stats

    def _getdata(self,index=slice(None)):
        """
        Args:
            index (slice or num):
        Returns:
            array: nlines x nspec x nchan x ndet
        """
        lines = self._getlines()
        self._needlines(lines)

        if index is not Ellipsis:
            lines = lines.__getitem__(index)
        return np.stack([line.data for line in lines],axis=0)

    def _getstats(self,index=slice(None)):
        """
        Args:
            index (slice or num):
        Returns:
            array: nlines x nspec x nstats x ndet
        """
        lines = self._getlines()
        self._needlines(lines)

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

