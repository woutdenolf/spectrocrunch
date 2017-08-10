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

from .. import edf

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
        #import cProfile
        #self.pr = cProfile.Profile()
        #self.pr.enable()

        self.dir = TempDirectory()

    def tearDown(self):
        self.dir.cleanup()

        #self.pr.disable()
        #self.pr.dump_stats("keep.cprof")

    def test_nameparsing(self):
        paths = ['/tmp/a1','/tmp/a2','/tmp/b1']
        radix = ['a','a','b']
        mapnums = [range(0,100),range(100,200),range(0,50)]
        linenums = [range(0,100),range(0,100),range(0,50)]
        labels = [['ctr1','ctr2','xia00','xia01','xiaS0','xiast']]*3
        
        p = zip(paths,radix,mapnums,linenums,labels)

        # Ordered list
        validlabel = lambda label,linenum: linenum==0 if 'xia' not in label else True
        filesordered = [os.path.join(path,"{}_{}_{:04d}_0000_{:04d}.edf".format(radix,label,mapnum,linenum)) for path,radix,mapnums,linenums,labels in p for mapnum in mapnums for linenum in linenums for label in labels if validlabel(label,linenum)]

        # Grouped
        filesgrouped = collections.OrderedDict()
        for path,radix,mapnums,linenums,labels in p:
            if radix not in filesgrouped:
                filesgrouped[radix] = collections.OrderedDict()
            for mapnum in mapnums:
                if mapnum not in filesgrouped[radix]:
                    filesgrouped[radix][mapnum] = collections.OrderedDict()
                for linenum in linenums:
                    filesgrouped[radix][mapnum][linenum] = [os.path.join(path,"{}_{}_{:04d}_0000_{:04d}.edf".format(radix,label,mapnum,linenum))  for label in labels if validlabel(label,linenum)]

        # Randomize list of files
        files = copy(filesordered)
        shuffle(files)

        # Check sorting
        files2 = sorted(files,key=xiaedf.xiasortkey)
        print files2
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

    def test_memmap(self):
        path = self.dir.path
        radix = "memmap"
        mapnum = 0
        linenum = 0
        line = xiaedf.xialine_number(path,radix,mapnum,linenum)

        data = np.random.rand(8,5,1)
        line.save(data,["00"])

        emap = edf.edfmemmap(line.datafilenames()[0])
        fmap = edf.edffabio(line.datafilenames()[0])

        np.testing.assert_array_equal(data[...,0],fmap.data)
        np.testing.assert_array_equal(fmap.data,emap.data)

    def _testdata(self,dataorg,data,stats,xiaobject,dshape,sshape,ishape):
        # data.shape:  ....,nchan,ndet
        # stats.shape: ....,nstat,ndet

        ndim = len(dshape)
        
        # Check data
        xiaobject.onlyicrocr(False)
        xiaobject.dtcor(False)
        xiaobject.skipdetectors([])

        self.assertEqual(xiaobject.dshape,dshape)
        self.assertEqual(xiaobject.sshape,sshape)

        self.assertEqual(data.shape,xiaobject.dshape)
        self.assertEqual(stats.shape,xiaobject.sshape)
        self.assertEqual(data.dtype,xiaobject.dtype)
        self.assertEqual(stats.dtype,xiaobject.stype)

        np.testing.assert_array_equal(data,xiaobject.data)
        np.testing.assert_array_equal(stats,xiaobject.stats)

        # Check DT correction
        xiaobject.onlyicrocr(True)
        
        xiaobject.dtcor(False)
        xiaobject.norm("")
        icrocr = xiaobject.stats
        icr = icrocr[...,xiaobject.indexicr,:].reshape(ishape)
        ocr = icrocr[...,xiaobject.indexocr,:].reshape(ishape)
        cor = np.asarray(icr,dtype=xiaobject.CORTYPE) / np.asarray(ocr,dtype=xiaobject.CORTYPE)
        np.testing.assert_array_equal(dataorg[...,0],xiaobject.data*cor)
        
        xiaobject.onlyicrocr(False)
        xiaobject.dtcor(True)
        xiaobject.norm("")
        np.testing.assert_array_equal(dataorg[...,0],xiaobject.data)
        
        # Test slicing
        xiaobject.skipdetectors([])

        indices = genindexing.genindexingn(dshape,advanced=True,eco=True,nmax=50)#,nmax=50)
        for dsum in [False,True]:
            xiaobject.detectorsum(dsum)
            for norm in [False,True]:
                xiaobject.norm("flux" if norm else None)
                for dtcor in [False,True]:
                    xiaobject.dtcor(dtcor)
                    for onlyicrocr in [False,True]:
                        xiaobject.onlyicrocr(onlyicrocr)
                        for together in [False,True]:
                            
                            for index in indices:
                                print "\n"*10
                                #index = (slice(8, -5, 4), slice(2, -8, 1), [-616, 854], 1)
                                #index = (Ellipsis, -2, None, [False, False, True, True], None)
                                #index = (None,slice(None),)#(None, None, None, Ellipsis, None)
                                #index = (9, 13, [0,0,0], 0, slice(None))
                                #index = slice(None)
                                #index = ([6, -8], [-8, -13], -4, slice(96, -89, 347), slice(None, None, None))
                                #index = ([0,1], [0,1],0,0,slice(None, None, None))
                                #index = (-14, -8, [0,1], [0,0])
                                #index = (0,0, [1,2], [1,-1], slice(None, None, None))
                                if dsum:
                                    index = indexing.replace(index,ndim,[-1],[slice(None)])
                                    
                                print "================",index,"================"
                                print dsum,norm,dtcor,onlyicrocr,together

                                if together:
                                    xiaobject.dataandstats()
                                    ldata,lstats = xiaobject[index]
                                else:
                                    xiaobject.onlydata()
                                    ldata = xiaobject[index]
                                    xiaobject.onlystats()
                                    lstats = xiaobject[index]

                                # Get data directly
                                if dtcor and norm:
                                    ldata2 = dataorg[...,2]
                                elif dtcor:
                                    ldata2 = dataorg[...,0]
                                elif norm:
                                    ldata2 = dataorg[...,1]
                                else:
                                    ldata2 = data
                                
                                ldata2 = ldata2[index]
                                if dsum:
                                    ldata2 = xiaobject.sumdata(ldata2,xiaobject._getaxis(-1))

                                # Get stats directly
                                if xiaobject.nstats==2:
                                    lstats2 = stats[...,[xiaobject.STICR,xiaobject.STOCR],:]
                                else:
                                    lstats2 = stats
                                lstats2 = lstats2[indexing.replacefull(index,ndim,[-2])]
                                
                                # Check data
                                np.testing.assert_allclose(ldata,ldata2)

                                # Check stats
                                np.testing.assert_array_equal(lstats,lstats2)

                                #return

        # Test slicing vs. skip detector
        xiaobject.onlyicrocr(False)
        xiaobject.dtcor(False)
        xiaobject.detectorsum(False)
        xiaobject.norm(None)
        
        ndet = dshape[-1]

        if ndet>2:
            xiaobject.skipdetectors([])
            ind = range(0,ndet,2)

            xiaobject.dataandstats()
            ldata,lstats = xiaobject[...,ind]
            np.testing.assert_array_equal(data[...,ind],ldata)
            np.testing.assert_array_equal(stats[...,ind],lstats)

            xiaobject.skipdetectors([0,ndet-1])
            np.testing.assert_array_equal(data[...,1:-1],xiaobject.data)
            np.testing.assert_array_equal(stats[...,1:-1],xiaobject.stats)

        if ndet>1:
            xiaobject.skipdetectors([])

            xiaobject.dataandstats()
            ldata,lstats = xiaobject[...,1:]
            np.testing.assert_array_equal(data[...,1:],ldata)
            np.testing.assert_array_equal(stats[...,1:],lstats)

            xiaobject.skipdetectors([0])
            self.assertEqual(ldata.shape,xiaobject.dshape)
            self.assertEqual(lstats.shape,xiaobject.sshape)

            np.testing.assert_array_equal(data[...,1:],xiaobject.data)
            np.testing.assert_array_equal(stats[...,1:],xiaobject.stats)

    def _test_line(self,path,radix,mapnum,linenum,ndet,nspec,nchan):
        # Generate some spectra + statistics
        dataorg,data,stats = xiagen.data(nspec,nchan,ndet)

        # Save data
        line = xiaedf.xialine_number(path,radix,mapnum,linenum)
        xialabels = ["xia{:02d}".format(i) for i in range(ndet)]
        line.save(data,xialabels,stats=stats)

        # Check saved files
        expectedfiles = ["{}_xia{:02}_{:04}_0000_{:04}.edf".format(radix,det,mapnum,linenum) for det in range(ndet)]+\
                        ["{}_xiast_{:04}_0000_{:04}.edf".format(radix,mapnum,linenum)]
        expectedfiles.sort()
        self.dir.compare(expectedfiles,path=path)

        # Check static files names
        files = xiaedf.xiasearch(path,radix=radix)
        line2 = xiaedf.xialine_files(files)
        line3 = xiaedf.xialine_number(path,radix,mapnum,linenum)
        self.assertEqual(files,sorted([os.path.join(path,f) for f in expectedfiles],key=xiaedf.xiasortkey))
        self.assertEqual(line.statfilename(),line2.statfilename())
        self.assertEqual(line.datafilenames(),line2.datafilenames())
        self.assertEqual(line.statfilename(),line3.statfilename())
        self.assertEqual(line.datafilenames(),line3.datafilenames())

        # Check data
        dshape = (nspec,nchan,ndet)
        sshape = (nspec,xiaedf.xiadata.NSTATS,ndet)
        ishape = (nspec,1,ndet)
        self._testdata(dataorg,data,stats,line,dshape,sshape,ishape)
        #import matplotlib.pyplot as plt
        #plt.plot(data[0,:,0])
        #plt.show()

    def _xiaconfigidcheck(self,xiaobject):
        if isinstance(xiaobject,xiaedf.xiacompound):
            ret = []
            for l in xiaobject._items:
                add = self._xiaconfigidcheck(l)
                if isinstance(add,list):
                    ret += add
                else:
                    ret.append(add)
            return ret
        else:
            return (id(xiaobject._xiaconfig),id(xiaobject._cache))

    def _test_image(self,path,radix,mapnum,ndet,ncol,nrow,nchan):
        # Generate some spectra + statistics
        flux = np.linspace(1,2,nrow*ncol)
        dataorg,data,stats = xiagen.data(nrow*ncol,nchan,ndet,flux=flux)
        dataorg = dataorg.reshape(nrow,ncol,nchan,ndet,3)
        data = data.reshape(nrow,ncol,nchan,ndet)
        stats = stats.reshape(nrow,ncol,xiaedf.xiadata.NSTATS,ndet)
        ctrs = {"flux":flux.reshape(nrow,ncol)}

        # Save data
        image = xiaedf.xiaimage_number(path,radix,mapnum)
        xialabels = ["xia{:02d}".format(i) for i in range(ndet)]
        image.save(data,xialabels,stats=stats,ctrs=ctrs)

        # Check saved files
        expectedfiles = ["{}_xia{:02}_{:04}_0000_{:04}.edf".format(radix,det,mapnum,linenum) for det in range(ndet) for linenum in range(nrow)]+\
                        ["{}_xiast_{:04}_0000_{:04}.edf".format(radix,mapnum,linenum) for linenum in range(nrow)]+\
                        ["{}_{}_{:04}_0000.edf".format(radix,k,mapnum) for k in ctrs]
                        
        expectedfiles.sort()
        self.dir.compare(expectedfiles,path=path)

        # Check files names
        files = xiaedf.xiasearch(path,radix=radix)
        image2 = xiaedf.xiaimage_files(files)
        image3 = xiaedf.xiaimage_linenumbers(path,radix,mapnum,range(nrow))
        image4 = xiaedf.xiaimage_number(path,radix,mapnum)
        self.assertEqual(files,sorted([os.path.join(path,f) for f in expectedfiles],key=xiaedf.xiasortkey))
        
        self.assertEqual(image.statfilenames(),image2.statfilenames())
        self.assertEqual(image.datafilenames(),image2.datafilenames())
        self.assertEqual(image.ctrfilenames(),image2.ctrfilenames())
        
        self.assertEqual(image.statfilenames(),image3.statfilenames())
        self.assertEqual(image.datafilenames(),image3.datafilenames())
        self.assertEqual(image.ctrfilenames(),image3.ctrfilenames())
        
        self.assertEqual(image.statfilenames(),image4.statfilenames())
        self.assertEqual(image.datafilenames(),image4.datafilenames())
        self.assertEqual(image.ctrfilenames(),image4.ctrfilenames())
        
        # Check only one config
        for o in [image,image2,image3,image4]:
            tmp = self._xiaconfigidcheck(o)
            self.assertEqual(tmp.count(tmp[0]),len(tmp))

        # Check data
        dshape = (nrow,ncol,nchan,ndet)
        sshape = (nrow,ncol,xiaedf.xiadata.NSTATS,ndet)
        ishape = (nrow,ncol,1,ndet)
        self._testdata(dataorg,data,stats,image,dshape,sshape,ishape)

    def _test_stack(self,path,radix,ndet,ncol,nrow,nenergy,nchan):
        # Generate some spectra + statistics
        flux = np.linspace(1,2*nenergy,nrow*ncol*nenergy)
        dataorg,data,stats = xiagen.data(nrow*ncol*nenergy,nchan,ndet,flux=flux)
        dataorg = dataorg.reshape(nenergy,nrow,ncol,nchan,ndet,3)
        data = data.reshape(nenergy,nrow,ncol,nchan,ndet)
        stats = stats.reshape(nenergy,nrow,ncol,xiaedf.xiadata.NSTATS,ndet)
        ctrs = {"flux":flux.reshape(nenergy,nrow,ncol)}
        
        # Save data
        stack = xiaedf.xiastack_radix(path,radix)
        xialabels = ["xia{:02d}".format(i) for i in range(ndet)]
        stack.save(data,xialabels,stats=stats,ctrs=ctrs)

        # Check saved files
        expectedfiles = ["{}_xia{:02}_{:04}_0000_{:04}.edf".format(radix,det,mapnum,linenum) for det in range(ndet) for linenum in range(nrow) for mapnum in range(nenergy)]+\
                        ["{}_xiast_{:04}_0000_{:04}.edf".format(radix,mapnum,linenum) for linenum in range(nrow) for mapnum in range(nenergy)]+\
                        ["{}_{}_{:04}_0000.edf".format(radix,k,mapnum) for mapnum in range(nenergy) for k in ctrs]
                        
        expectedfiles.sort()
        self.dir.compare(expectedfiles,path=path)

        # Check files names
        files = xiaedf.xiasearch(path,radix=radix)
        stack2 = xiaedf.xiastack_files(files)
        stack3 = xiaedf.xiastack_radix(path,radix)
        self.assertEqual(files,sorted([os.path.join(path,f) for f in expectedfiles],key=xiaedf.xiasortkey))
        self.assertEqual(stack.statfilenames(),stack2.statfilenames())
        self.assertEqual(stack.datafilenames(),stack2.datafilenames())
        self.assertEqual(stack.ctrfilenames(),stack2.ctrfilenames())
        
        self.assertEqual(stack.statfilenames(),stack3.statfilenames())
        self.assertEqual(stack.ctrfilenames(),stack3.ctrfilenames())

        # Check data
        dshape = (nenergy,nrow,ncol,nchan,ndet)
        sshape = (nenergy,nrow,ncol,xiaedf.xiadata.NSTATS,ndet)
        ishape = (nenergy,nrow,ncol,1,ndet)
        self._testdata(dataorg,data,stats,stack,dshape,sshape,ishape)

    def test_line(self):
        mapnum = 2
        linenum = 10

        i = 0
        for ndet in [4,1]:
            for nchan in [1024,1]:
                for nspec in [10,1]:
                    self._test_line(os.path.join(self.dir.path,"test_line_{}".format(i)),"test_line_{}".format(i),mapnum,linenum,ndet,nspec,nchan)
                    i += 1

    def test_image(self):
        mapnum = 2

        i = 0
        for ndet in [4,1]:
            for nchan in [1024,1]:
                for ncol in [10,1]:
                    for nrow in [15,1]:
                        self._test_image(os.path.join(self.dir.path,"test_image_{}".format(i)),"test_image_{}".format(i),mapnum,ndet,ncol,nrow,nchan)
                        i += 1

    def test_stack(self):

        i = 0
        for ndet in [4]:
            for nchan in [1024]:
                for ncol in [10]:
                    for nrow in [15]:
                        for nenergy in [20]:
                            self._test_stack(os.path.join(self.dir.path,"test_stack_{}".format(i)),"test_stack_{}".format(i),ndet,ncol,nrow,nenergy,nchan)
                            i += 1



def test_suite_all():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    #testSuite.addTest(test_xiaedf("test_memmap"))
    #testSuite.addTest(test_xiaedf("test_nameparsing"))
    #testSuite.addTest(test_xiaedf("test_line"))
    testSuite.addTest(test_xiaedf("test_image"))
    #testSuite.addTest(test_xiaedf("test_stack"))
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite_all()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
