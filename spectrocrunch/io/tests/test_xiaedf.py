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

import unittest

from .. import xiaedf

from . import xiagen

from testfixtures import TempDirectory

import numpy as np

import os

from glob import glob

from random import shuffle

import collections

from copy import copy

from ...common import indexing

from ...common.tests import genindexing





class test_xiaedf(unittest.TestCase):

    def setUp(self):
        import cProfile
        self.pr = cProfile.Profile()
        self.pr.enable()

        self.dir = TempDirectory()

    def tearDown(self):
        self.dir.cleanup()

        self.pr.disable()
        self.pr.dump_stats("keep.cprof")

    def test_nameparsing(self):
        paths = ['/tmp/a1','/tmp/a2','/tmp/b1']
        radix = ['a','a','b']
        mapnums = [range(0,100),range(100,200),range(0,50)]
        linenums = [range(0,100),range(0,100),range(0,50)]
        labels = [['00','01','S0','st']]*3
        
        p = zip(paths,radix,mapnums,linenums,labels)

        # Ordered list
        filesordered = [os.path.join(path,"{}_xia{}_{:04d}_0000_{:04d}.edf".format(radix,label,mapnum,linenum))  for path,radix,mapnums,linenums,labels in p for mapnum in mapnums for linenum in linenums for label in labels]

        # Grouped
        filesgrouped = collections.OrderedDict()
        for path,radix,mapnums,linenums,labels in p:
            if radix not in filesgrouped:
                filesgrouped[radix] = collections.OrderedDict()
            for mapnum in mapnums:
                if mapnum not in filesgrouped[radix]:
                    filesgrouped[radix][mapnum] = collections.OrderedDict()
                for linenum in linenums:
                    filesgrouped[radix][mapnum][linenum] = [os.path.join(path,"{}_xia{}_{:04d}_0000_{:04d}.edf".format(radix,label,mapnum,linenum))  for label in labels]

        # Randomize list of files
        files = copy(filesordered)
        shuffle(files)

        # Check sorting
        files2 = sorted(files,key=xiaedf.xiasortkey)
        self.assertEqual(filesordered,files2)

        # Check grouping
        g = xiaedf.xiagroup(files)
        self.assertEqual(filesgrouped,g)

        # Check map grouping
        shuffle(files)
        files = [f for f in files if os.path.basename(f).startswith('a')]
        g = xiaedf.xiagroupmaps(files)
        self.assertEqual(filesgrouped['a'],g)

        # Check line grouping
        shuffle(files)
        files = [f for f in files if '0100_0000' in os.path.basename(f)]
        g = xiaedf.xiagrouplines(files)
        self.assertEqual(filesgrouped['a'][100],g)

    def _testdata(self,dataorg,data,stats,o,dshape,sshape,ishape):
        # data.shape:  ....,nchan,ndet
        # stats.shape: ....,nstat,ndet

        ndim = len(dshape)

        # Check data
        o.onlyicrocr(False)
        o.dtcor(False)
        o.skipdetectors([])

        self.assertEqual(o.dshape,dshape)
        self.assertEqual(o.sshape,sshape)

        self.assertEqual(data.shape,o.dshape)
        self.assertEqual(stats.shape,o.sshape)
        self.assertEqual(data.dtype,o.dtype)
        self.assertEqual(stats.dtype,o.stype)

        np.testing.assert_array_equal(data,o.data)
        np.testing.assert_array_equal(stats,o.stats)

        # Check DT correction
        o.onlyicrocr(True)
        icrocr = o.stats
        np.testing.assert_allclose(dataorg,o.data*icrocr[...,0,:].reshape(ishape)/icrocr[...,1,:].reshape(ishape))
        o.onlyicrocr(False)
        o.dtcor(True)
        np.testing.assert_allclose(dataorg,o.data)
        
        # Test slicing
        o.skipdetectors([])
        o.dataandstats()
        o.onlyicrocr(False)
        o.dtcor(True)

        indices = genindexing.genindexingn(dshape,advanced=True)
        for index in indices:
            ldata,lstats = o[index]
            ldata2 = dataorg[index]

            np.testing.assert_array_equal(ldata,ldata2)

            index2 = indexing.replacefull(index,ndim,-2)
            lstats2 = stats[index2]

            np.testing.assert_array_equal(lstats,lstats2)

        # Test slicing vs. skip detector
        o.dataandstats()

        ndet = dshape[-1]
        if ndet>1:
            o.skipdetectors([])
            ldata,lstats = o[...,1:]
            np.testing.assert_array_equal(data[...,1:],ldata)
            np.testing.assert_array_equal(stats[...,1:],lstats)

            o.skipdetectors([0])
            self.assertEqual(ldata.shape,o.dshape)
            self.assertEqual(lstats.shape,o.sshape)

            np.testing.assert_array_equal(data[...,1:],o.data)
            np.testing.assert_array_equal(stats[...,1:],o.stats)

        if ndet>2:
            o.skipdetectors([])
            ind = range(0,ndet,2)

            ldata,lstats = o[...,ind]
            np.testing.assert_array_equal(data[...,ind],ldata)
            np.testing.assert_array_equal(stats[...,ind],lstats)

            o.skipdetectors([0,ndet-1])
            np.testing.assert_array_equal(data[...,1:-1],o.data)
            np.testing.assert_array_equal(stats[...,1:-1],o.stats)

    def _test_line(self,path,radix,mapnum,linenum,ndet,nspec,nchan):
        # Generate some spectra + statistics
        dataorg,data,stats = xiagen.data(nspec,nchan,ndet)

        # Save data
        line = xiaedf.xialine_number(path,radix,mapnum,linenum)
        xialabels = ["{:02d}".format(i) for i in range(ndet)]
        line.save(data,xialabels,stats=stats)

        # Check saved files
        expectedfiles = ["{}_xia{:02}_{:04}_0000_{:04}.edf".format(radix,det,mapnum,linenum) for det in range(ndet)]+\
                        ["{}_xiast_{:04}_0000_{:04}.edf".format(radix,mapnum,linenum)]
        self.dir.compare(expectedfiles,path=path)

        # Check static files names
        files = xiaedf.xiasearch(path,radix=radix)
        line2 = xiaedf.xialine_files(files)
        self.assertEqual(files,sorted([os.path.join(path,f) for f in expectedfiles],key=xiaedf.xiasortkey))
        self.assertEqual(line.statfilename(),line2.statfilename())
        self.assertEqual(line.datafilenames(),line2.datafilenames())

        # Check data
        dshape = (nspec,nchan,ndet)
        sshape = (nspec,xiaedf.xiadata.NSTATS,ndet)
        ishape = (nspec,1,ndet)
        self._testdata(dataorg,data,stats,line,dshape,sshape,ishape)
        self._testdata(dataorg,data,stats,line2,dshape,sshape,ishape)
        #import matplotlib.pyplot as plt
        #plt.plot(data[0,:,0])
        #plt.show()

    def test_line(self):
        mapnum = 2
        linenum = 10

        i = 0
        for ndet in [1,4]:
            for nchan in [1,1024]:
                for nspec in [1,10]:
                    self._test_line(os.path.join(self.dir.path,"test_line_{}".format(i)),"test_line_{}".format(i),mapnum,linenum,ndet,nspec,nchan)
                    i += 1

    def _test_image(self,path,radix,mapnum,ndet,ncol,nrow,nchan):
        # Generate some spectra + statistics
        dataorg,data,stats = xiagen.data(nrow*ncol,nchan,ndet)
        dataorg = dataorg.reshape(nrow,ncol,nchan,ndet)
        data = data.reshape(nrow,ncol,nchan,ndet)
        stats = stats.reshape(nrow,ncol,xiaedf.xiadata.NSTATS,ndet)

        # Save data
        image = xiaedf.xiaimage_number(path,radix,mapnum)
        xialabels = ["{:02d}".format(i) for i in range(ndet)]
        image.save(data,xialabels,stats=stats)

        # Check saved files
        expectedfiles = ["{}_xia{:02}_{:04}_0000_{:04}.edf".format(radix,det,mapnum,linenum) for det in range(ndet) for linenum in range(nrow)]+\
                        ["{}_xiast_{:04}_0000_{:04}.edf".format(radix,mapnum,linenum) for linenum in range(nrow)]
        self.dir.compare(expectedfiles,path=path)

        # Check files names
        files = xiaedf.xiasearch(path,radix=radix)
        image2 = xiaedf.xiaimage_files(files)
        image3 = xiaedf.xiaimage_linenumbers(path,radix,mapnum,range(nrow))
        self.assertEqual(files,sorted([os.path.join(path,f) for f in expectedfiles],key=xiaedf.xiasortkey))
        self.assertEqual(image.statfilenames(),image2.statfilenames())
        self.assertEqual(image.datafilenames(),image2.datafilenames())
        self.assertEqual(image.statfilenames(),image3.statfilenames())
        self.assertEqual(image.datafilenames(),image3.datafilenames())

        # Check data
        dshape = (nrow,ncol,nchan,ndet)
        sshape = (nrow,ncol,xiaedf.xiadata.NSTATS,ndet)
        ishape = (nrow,ncol,1,ndet)
        self._testdata(dataorg,data,stats,image,dshape,sshape,ishape)
        self._testdata(dataorg,data,stats,image2,dshape,sshape,ishape)

    def test_image(self):
        mapnum = 2

        i = 0
        for ndet in [1,4]:
            for nchan in [1,1024]:
                for ncol in [1,10]:
                    for nrow in [1,10]:
                        self._test_image(os.path.join(self.dir.path,"test_image_{}".format(i)),"test_image_{}".format(i),mapnum,ndet,ncol,nrow,nchan)
                        i += 1


    def test_memmap(self):
        path = self.dir.path
        radix = "memmap"
        mapnum = 0
        linenum = 0
        line = xiaedf.xialine_number(path,radix,mapnum,linenum)

        data = np.random.rand(8,5,1)
        line.save(data,["00"])

        emap = xiaedf.edfmemmap(line.datafilenames()[0])

        np.testing.assert_array_equal(data[...,0],emap.data)

def test_suite_all():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    #testSuite.addTest(test_xiaedf("test_memmap"))
    #testSuite.addTest(test_xiaedf("test_nameparsing"))
    testSuite.addTest(test_xiaedf("test_line"))
    #testSuite.addTest(test_xiaedf("test_image"))
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite_all()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
