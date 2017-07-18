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

import logging
logger = logging.getLogger(__name__)

class xialine(object):

    STDET = 0
    STEVT = 1
    STICR = 2
    STOCR = 3
    STLT = 4
    STDT = 5
    NSTATS = 6
    
    def __init__(self,path,radix,mapnum,linenum,overwrite=False,skiponread=None):
        """
        Args:
            path(list): path
            radix(list): radix
            num(numbers.Integral): map number
            linenum(numbers.Integral): line number
            linenum(Optional(overwrite)): line number
            skiponread(Optional(list)): detector numbers to skip when reading
        Returns:
            None
        """
        
        self.filename = os.path.join(path,"{}_xia{{}}_{:04d}_0000_{:04d}.edf".format(radix,mapnum,linenum))
        self.overwrite = overwrite
        self.skiponread(skiponread)

    @property
    def data(self):
        return self._getdata()

    @property
    def stats(self):
        return self._getstats()

    @property
    def icrocr(self):
        stats = self._getstats()
        return stats[xialine.STICR,...],stats[xialine.STOCR,...]

    def __getitem__(self, index):
        data = self._getdata()
        stats = self._getstats()

        data = data.__getitem__(index)
        if stats is not None:
            stats = stats.__getitem__(index)

        return data,stats

    def skiponread(self,dets):
        """
        Args:
            dets(list): detector to skip when reading
        Returns:
            list: formated labels
        """
        if dets is None:
            self._skiponread = []
        else:
            self._skiponread = dets

    def _skiplabels(self):
        """
        Args:
            None
        Returns:
            list: formated labels
        """
        labels = copy(self._skiponread)
        for i in range(len(labels)):
            if isinstance(labels[i],numbers.Integral):
                labels[i] = "{:02}".format(labels[i])
        return labels

    def _skipnums(self):
        """
        Args:
            None
        Returns:
            list: detector numbers
        """

        nums = []
        for i in range(len(self._skiponread)):
            if isinstance(self._skiponread[i],str):
                if self._skiponread[i].isnumeric(): #"00", "01", ... (there are no statistics for sums)
                    nums.append(int(self._skiponread[i]))
            else:
                nums.append(self._skiponread[i])
        return nums
 
    def _getdata(self):
        """
        Args:
            None
        Returns:
            array: nspec x nchan x ndet
        """

        # Files and their detector labels
        files = glob(self.filename.format("*[0-9]*"))
        if len(files)==0:
            return None
        labels = [os.path.basename(f).split('_')[-4][3:] for f in files]

        # Sort on labels
        lst = zip(labels,files)
        lst.sort()

        # Skip when required
        skiplabels = self._skiplabels()
        files = [f for det,f in lst if det not in skiplabels]

        data = np.dstack([fabio.open(f).data for f in files])
        
        return data
    
    def _getstats(self):
        """
        Args:
            None
        Returns:
            array: nspec x nstats x ndet
        """
        f = self.filename.format("st")
        if os.path.isfile(f):
            stats = fabio.open(f).data
            s = stats.shape
            nspec = s[0]
            ndet = s[1]//xialine.NSTATS
            stats = stats.reshape((nspec,xialine.NSTATS,ndet))

            skipnums = self._skipnums()
            if len(skipnums)!=0:
                ind = [i for i in range(ndet) if i not in skipnums]
                stats = stats[...,ind]
            
        else:
            stats = None

        return stats

    def _write(self,img,filename):
        if not self.overwrite:
            if os.path.isfile(filename):
                logger.warning("{} not saved (already exists)".format(filename))
                return

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
            self._write(img,self.filename.format(xialabels[i]))

        if stats is not None:
            img = fabio.edfimage.EdfImage(data=stats.reshape(( stats.shape[0], stats.shape[1] * stats.shape[2] )))
            self._write(img,self.filename.format("st"))

class xiamap(object):

    def __init__(self,path,radix,num):
        """
        Args:
            path(list): path
            radix(list): radix
            num(numbers.Integral): map number

        Returns:
            None
        """
        self.filemask = os.path.join(path,"%s_xia*_%04d_0000_*.edf"%(radix,num))

class xiaedfnd(object):

    def __init__(self,datadirs,scannames,scannumbers):
        """
        Args:
            datadirs(list): path
            scannames(list): radix
            scannumbers(list(list(numbers.Integral))): scan numbers

        Returns:
            None
        """

        if isinstance(scannumbers,numbers.Integral):
            scannumbers = [[scannumbers]]
        elif isinstance(scannumbers,list):
            if isinstance(scannumbers[0],numbers.Integral):
                scannumbers = [scannumbers]

        n = len(scannumbers)
        if isinstance(datadirs,str):
            datadirs = [datadirs]*n
        if isinstance(scannames,str):
            scannames = [scannames]*n

        #for path,name,num:
        #    # [path]_[name]_xia[detnum]_[nume]_0000_[linenumber].edf
        #    filemask = os.path.join(path,"%s_xia*_%04d_0000_*.edf"%(name,num))

        
