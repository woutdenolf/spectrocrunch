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
from testfixtures import TempDirectory
import numpy as np
import os
from glob import glob
from random import shuffle
import collections
from copy import copy

from .. import xiaedf
from .. import edf
from . import xiagen
from ...utils import indexing
from ...utils import listtools
from ...utils.tests import genindexing

import logging
logger = logging.getLogger(__name__)

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

    def test_nameparsing_special(self):
        self.assertEqual(xiaedf.xianameparser.parse("test_puz_PUZ_xmap_x1_00_0009_0000.edf"),\
        xiaedf.XiaName(radix="test_puz",mapnum=9,linenum=-1,label="PUZ_xmap_x1_00",baselabel="PUZ_xmap_x1",detector="00"))
 
        self.assertEqual(xiaedf.xianameparser.parse("samB6_mapa_xmap_x3c_00_0002_0000.edf"),\
        xiaedf.XiaName(radix="samB6_mapa",mapnum=2,linenum=-1,label="xmap_x3c_00",baselabel="xmap_x3c",detector="00"))
        
        self.assertEqual(xiaedf.xianameparser.parse("samB6_mapa_xiast_0002_0000_0069.edf"),\
        xiaedf.XiaName(radix="samB6_mapa",mapnum=2,linenum=69,label="xiast",baselabel="xia",detector="st"))

    def test_nameparsing(self):
        paths = ['/tmp/a1','/tmp/a2','/tmp/b1']
        radix = ['a','a','b']
        mapnums = [range(0,100),range(100,200),range(0,50)]
        linenums = [range(0,100),range(0,100),range(0,50)]
        labels = [['arr_1','arr_2','xia00','xia01','xiaS0','xiast']]*3
        
        p = zip(paths,radix,mapnums,linenums,labels)

        # Grouped
        filesgrouped = collections.OrderedDict()
        for path,radix,mapnums,linenums,labels in p:
            if radix not in filesgrouped:
                filesgrouped[radix] = collections.OrderedDict()
            for mapnum in mapnums:
                if mapnum not in filesgrouped[radix]:
                    filesgrouped[radix][mapnum] = collections.OrderedDict()
                linenum = -1
                filesgrouped[radix][mapnum][linenum] = [os.path.join(path,xiaedf.xiafilename(radix,mapnum,linenum,label)) for label in labels if "xia" not in label]
                for linenum in linenums:
                    filesgrouped[radix][mapnum][linenum] = [os.path.join(path,xiaedf.xiafilename(radix,mapnum,linenum,label)) for label in labels if "xia" in label]
 
        # Ordered list
        filesordered = [vline for radix,vradix in filesgrouped.items() for mapnum,vmap in vradix.items() for linenum,vline in vmap.items()]
        filesordered = list(listtools.flatten(filesordered))
 
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

    def test_memmap(self):
        path = self.dir.path
        radix = "memmap"
        mapnum = 0
        linenum = 0
        line = xiaedf.xialine_number(path,radix,mapnum,linenum)

        data = np.random.rand(8,5,1)
        line.save(data,["xia00"])

        emap = edf.edfmemmap(line.datafilenames()[0])
        fmap = edf.edfimage(line.datafilenames()[0])

        np.testing.assert_array_equal(data[...,0],fmap.data)
        np.testing.assert_array_equal(fmap.data,emap.data)

    def _testdata(self,dataorg,data,stats,ctrs,xiaobject,dshape,sshape,ishape,cshape):
        # data.shape:  ....,nchan,ndet
        # stats.shape: ....,nstat,ndet

        ndim = len(dshape)
        
        # Check data
        xiaobject.onlyicrocr(False)
        xiaobject.dtcor(False)
        xiaobject.detectorsum(False)
        xiaobject.globalnorm(None)
        xiaobject.exclude_detectors([])

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
        icrocr = xiaobject.stats
        icr = icrocr[...,xiaobject.indexicr,:].reshape(ishape)
        ocr = icrocr[...,xiaobject.indexocr,:].reshape(ishape)
        cor = np.asarray(icr,dtype=xiaobject.CORTYPE) / np.asarray(ocr,dtype=xiaobject.CORTYPE)
        np.testing.assert_array_equal(dataorg[...,0],xiaobject.data*cor)
        
        xiaobject.onlyicrocr(False)
        xiaobject.dtcor(True)
        np.testing.assert_array_equal(dataorg[...,0],xiaobject.data)

        # Test slicing
        xiaobject.exclude_detectors([])
 
        indices = genindexing.genindexingn(dshape,advanced=True,eco=True,nmax=50)
        for dsum in [False,True]:
            xiaobject.detectorsum(dsum)
            for norm in [False,True]:
                xiaobject.globalnorm("arr_flux" if norm else None)
                for dtcor in [False,True]:
                    xiaobject.dtcor(dtcor)
                    for onlyicrocr in [False,True]:
                        xiaobject.onlyicrocr(onlyicrocr)
                        for together in [False,True]:
                            
                            for index in indices:
                                #index = (-10, Ellipsis, slice(8, -3, 2), [True, True, False, False])

                                if dsum:
                                    index = indexing.replace(index,ndim,[-1],[slice(None)])

                                logger.debug("\n"*5)
                                logger.debug("index = {}".format(index))
                                logger.debug("sum = {}".format(dsum))
                                logger.debug("norm = {}".format(norm))
                                logger.debug("dtcor = {}".format(dtcor))
                                logger.debug("onlyicrocr = {}".format(onlyicrocr))
                                logger.debug("together = {}".format(together))

                                if together:
                                    if ctrs is None:
                                        xiaobject.dataandstats()
                                        ldata,lstats = xiaobject[index]
                                    else:
                                        xiaobject.dataall()
                                        ldata,lstats,lcounters = xiaobject[index]
                                else:
                                    xiaobject.onlydata()
                                    ldata = xiaobject[index]
                                    xiaobject.onlystats()
                                    lstats = xiaobject[index]
                                    if ctrs is not None:
                                        xiaobject.onlycounters()
                                        lcounters = xiaobject[index]
                                
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

                                # Check data
                                np.testing.assert_allclose(ldata,ldata2)
                                
                                # Get stats directly
                                if xiaobject.nstats==2:
                                    lstats2 = stats[...,[xiaobject.STICR,xiaobject.STOCR],:]
                                else:
                                    lstats2 = stats
                                lstats2 = lstats2[indexing.replacefull(index,ndim,[-2])]
                                
                                # Check stats
                                np.testing.assert_array_equal(lstats,lstats2)
                                
                                # Get counters directly
                                if ctrs is not None:
                                    try:
                                        lcounters2 = ctrs[indexing.replacefull(index,ndim,[-2,-1])]
                                        check = True
                                    except:
                                        # This happens when the MCA channel index cannot be applied to the counter index (same dimension)
                                        check = False
                                    
                                    # Check counters
                                    if check:
                                        np.testing.assert_array_equal(lcounters,lcounters2)
                                
        # Test slicing vs. skip detector
        xiaobject.onlyicrocr(False)
        xiaobject.dtcor(False)
        xiaobject.detectorsum(False)
        xiaobject.globalnorm(None)
        
        ndet = dshape[-1]

        if ndet>2:
            xiaobject.exclude_detectors([])
            ind = range(0,ndet,2)

            xiaobject.dataandstats()
            ldata,lstats = xiaobject[...,ind]
            np.testing.assert_array_equal(data[...,ind],ldata)
            np.testing.assert_array_equal(stats[...,ind],lstats)

            xiaobject.exclude_detectors([0,ndet-1])
            np.testing.assert_array_equal(data[...,1:-1],xiaobject.data)
            np.testing.assert_array_equal(stats[...,1:-1],xiaobject.stats)

        if ndet>1:
            xiaobject.exclude_detectors([])

            xiaobject.dataandstats()
            ldata,lstats = xiaobject[...,1:]
            np.testing.assert_array_equal(data[...,1:],ldata)
            np.testing.assert_array_equal(stats[...,1:],lstats)

            xiaobject.exclude_detectors([0])
            self.assertEqual(ldata.shape,xiaobject.dshape)
            self.assertEqual(lstats.shape,xiaobject.sshape)

            np.testing.assert_array_equal(data[...,1:],xiaobject.data)
            np.testing.assert_array_equal(stats[...,1:],xiaobject.stats)

    def _testcopy(self,xiaobject,xiaobjectgen,path):
        iscompound = isinstance(xiaobject,xiaedf.xiacompound)
    
        i = 0
        for dsum in [False,True]:
            xiaobject.detectorsum(dsum)
            copystats = not dsum
            if dsum:
                xialabels = ["xiaS1"]
            else:
                xialabels = xiaobject.xialabels_used
            for norm in [False,True]:
                xiaobject.globalnorm("arr_flux" if norm else None)
                for dtcor in [False,True]:
                    xiaobject.dtcor(dtcor)
                    for onlyicrocr in [False,True]:
                        xiaobject.onlyicrocr(onlyicrocr)
                        for copy in [False,True]:
                            for copyctrs in [True,False]:
                                for deflabels in [False,True]:
                                    if copyctrs and not iscompound: # no counter when not a xiacompound (counters are saved per image)
                                        continue
                                    path2 = os.path.join(path,"copy{}".format(i))
                                    i += 1
                                    xiaobject2 = xiaobjectgen(path2)
                                    
                                    logger.debug("\n"*5)
                                    logger.debug("sum = {}".format(dsum))
                                    logger.debug("norm = {}".format(norm))
                                    logger.debug("dtcor = {}".format(dtcor))
                                    logger.debug("onlyicrocr = {}".format(onlyicrocr))
                                    logger.debug("copyctrs = {}".format(copyctrs))
                                    logger.debug("deflabels = {}".format(deflabels))
                                    logger.debug("copy = {}".format(copy))

                                    kwargs={}
                                    if deflabels:
                                        kwargs["xialabels"] = xialabels
                                    if copy:
                                        kwargs["stats"] = copystats
                                        if iscompound:
                                            kwargs["ctrs"] = copyctrs
                                            if copyctrs:
                                                kwargs["ctrnames"] = xiaobject.counterbasenames()
                                        xiaobject2.save(xiaobject,**kwargs)
                                    else:
                                        if copystats:
                                            kwargs["stats"] = xiaobject.allstats
                                        if copyctrs and iscompound:
                                            kwargs["ctrs"] = xiaobject.counters[...,0]
                                            kwargs["ctrnames"] = xiaobject.counterbasenames()
                                        xiaobject2.save(xiaobject.data,**kwargs)

                                    np.testing.assert_array_equal(xiaobject.data,xiaobject2.data)
                                    if copystats:
                                        np.testing.assert_array_equal(xiaobject.allstats,xiaobject2.allstats)
                                    if copyctrs:
                                        np.testing.assert_array_equal(xiaobject.counters,xiaobject2.counters)

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

        # Check copy
        linegen = lambda x: xiaedf.xialine_number(x,radix,mapnum,linenum)
        self._testcopy(line,linegen,path)

        # Check data
        dshape = (nspec,nchan,ndet)
        sshape = (nspec,xiaedf.xiadata.NSTATS,ndet)
        ishape = (nspec,1,ndet)
        self._testdata(dataorg,data,stats,None,line,dshape,sshape,ishape,None)
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
        ctr1 = np.linspace(10,20,nrow*ncol)
        ctr2 = np.linspace(20,30,nrow*ncol)
        dataorg,data,stats = xiagen.data(nrow*ncol,nchan,ndet,flux=flux)
        dataorg = dataorg.reshape(nrow,ncol,nchan,ndet,3)
        data = data.reshape(nrow,ncol,nchan,ndet)
        stats = stats.reshape(nrow,ncol,xiaedf.xiadata.NSTATS,ndet)
        ctrs = np.stack([ctr1,flux,ctr2],axis=1).reshape(nrow,ncol,3)
        ctrnames = ["arr_ctr1","arr_flux","arr_ctr2"]
        
        # Save data
        image = xiaedf.xiaimage_number(path,radix,mapnum)
        xialabels = ["xia{:02d}".format(i) for i in range(ndet)]
        image.save(data,xialabels,stats=stats,ctrs=ctrs,ctrnames=ctrnames)

        # Check saved files
        expectedfiles = ["{}_xia{:02}_{:04}_0000_{:04}.edf".format(radix,det,mapnum,linenum) for det in range(ndet) for linenum in range(nrow)]+\
                        ["{}_xiast_{:04}_0000_{:04}.edf".format(radix,mapnum,linenum) for linenum in range(nrow)]+\
                        ["{}_{}_{:04}_0000.edf".format(radix,k,mapnum) for k in ctrnames]
                        
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

        # Check copy
        imagegen = lambda x: xiaedf.xiaimage_number(x,radix,mapnum)
        self._testcopy(image,imagegen,path)
        
        # Check data
        dshape = (nrow,ncol,nchan,ndet)
        sshape = (nrow,ncol,xiaedf.xiadata.NSTATS,ndet)
        ishape = (nrow,ncol,1,ndet)
        cshape = (nrow,ncol,3,1)
        ctrs = ctrs[:,:,np.argsort(ctrnames)]
        ctrs = ctrs[...,np.newaxis]
        self._testdata(dataorg,data,stats,ctrs,image,dshape,sshape,ishape,cshape)

    def _test_stack(self,path,radix,ndet,ncol,nrow,nenergy,nchan):
        # Generate some spectra + statistics
        flux = np.linspace(1,2*nenergy,nrow*ncol*nenergy)
        ctr1 = np.linspace(10,20,nrow*ncol*nenergy)
        ctr2 = np.linspace(20,30,nrow*ncol*nenergy)
        dataorg,data,stats = xiagen.data(nrow*ncol*nenergy,nchan,ndet,flux=flux)
        dataorg = dataorg.reshape(nenergy,nrow,ncol,nchan,ndet,3)
        data = data.reshape(nenergy,nrow,ncol,nchan,ndet)
        stats = stats.reshape(nenergy,nrow,ncol,xiaedf.xiadata.NSTATS,ndet)
        ctrs = np.stack([ctr1,flux,ctr2],axis=1).reshape(nenergy,nrow,ncol,3)
        ctrnames = ["arr_ctr1","arr_flux","arr_ctr2"]
        
        # Save data
        stack = xiaedf.xiastack_radix(path,radix)
        xialabels = ["xia{:02d}".format(i) for i in range(ndet)]
        stack.save(data,xialabels,stats=stats,ctrs=ctrs,ctrnames=ctrnames)

        # Check saved files
        expectedfiles = ["{}_xia{:02}_{:04}_0000_{:04}.edf".format(radix,det,mapnum,linenum) for det in range(ndet) for linenum in range(nrow) for mapnum in range(nenergy)]+\
                        ["{}_xiast_{:04}_0000_{:04}.edf".format(radix,mapnum,linenum) for linenum in range(nrow) for mapnum in range(nenergy)]+\
                        ["{}_{}_{:04}_0000.edf".format(radix,k,mapnum) for mapnum in range(nenergy) for k in ctrnames]
                        
        expectedfiles.sort()
        self.dir.compare(expectedfiles,path=path)

        # Check files names
        files = xiaedf.xiasearch(path,radix=radix)
        stack2 = xiaedf.xiastack_files(files)
        stack3 = xiaedf.xiastack_radix(path,radix)
        stack4 = xiaedf.xiastack_mapnumbers(path,radix,range(nenergy))
        
        self.assertEqual(files,sorted([os.path.join(path,f) for f in expectedfiles],key=xiaedf.xiasortkey))
        self.assertEqual(stack.statfilenames(),stack2.statfilenames())
        self.assertEqual(stack.datafilenames(),stack2.datafilenames())
        self.assertEqual(stack.ctrfilenames(),stack2.ctrfilenames())
        
        self.assertEqual(stack.statfilenames(),stack3.statfilenames())
        self.assertEqual(stack.datafilenames(),stack3.datafilenames())
        self.assertEqual(stack.ctrfilenames(),stack3.ctrfilenames())

        self.assertEqual(stack.statfilenames(),stack4.statfilenames())
        self.assertEqual(stack.datafilenames(),stack4.datafilenames())
        self.assertEqual(stack.ctrfilenames(),stack4.ctrfilenames())

        # Check copy
        stackgen = lambda x: xiaedf.xiastack_radix(x,radix)
        self._testcopy(stack,stackgen,path)

        # Check data
        dshape = (nenergy,nrow,ncol,nchan,ndet)
        sshape = (nenergy,nrow,ncol,xiaedf.xiadata.NSTATS,ndet)
        ishape = (nenergy,nrow,ncol,1,ndet)
        cshape = (nenergy,nrow,ncol,3,1)
        ctrs = ctrs[:,:,:,np.argsort(ctrnames)]
        ctrs = ctrs[...,np.newaxis]
        self._testdata(dataorg,data,stats,ctrs,stack,dshape,sshape,ishape,cshape)

    def test_line(self):
        mapnum = 2
        linenum = 10

        i = 0
        for ndet in [4,1]:
            for nchan in [32,1]:
                for nspec in [10,1]:
                    self._test_line(os.path.join(self.dir.path,"test_line_{}".format(i)),"test_line_{}".format(i),mapnum,linenum,ndet,nspec,nchan)
                    i += 1

    def test_image(self):
        mapnum = 2

        i = 0
        for ndet in [4,1]:
            for nchan in [32,1]:
                for ncol in [10,1]:
                    for nrow in [15,1]:
                        self._test_image(os.path.join(self.dir.path,"test_image_{}".format(i)),"test_image_{}".format(i),mapnum,ndet,ncol,nrow,nchan)
                        i += 1

    def test_stack(self):

        i = 0
        for ndet in [4]:
            for nchan in [32]:
                for ncol in [10]:
                    for nrow in [15]:
                        for nenergy in [8]:
                            self._test_stack(os.path.join(self.dir.path,"test_stack_{}".format(i)),"test_stack_{}".format(i),ndet,ncol,nrow,nenergy,nchan)
                            i += 1


def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_xiaedf("test_memmap"))
    testSuite.addTest(test_xiaedf("test_nameparsing"))
    testSuite.addTest(test_xiaedf("test_nameparsing_special"))
    testSuite.addTest(test_xiaedf("test_line"))
    testSuite.addTest(test_xiaedf("test_image"))
    testSuite.addTest(test_xiaedf("test_stack"))
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
