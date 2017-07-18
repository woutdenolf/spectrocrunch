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

from .. import xiaedfnd

from testfixtures import TempDirectory

import numpy as np

class test_xiaedfnd(unittest.TestCase):

    def setUp(self):
        self.dir = TempDirectory()

    def tearDown(self):
        self.dir.cleanup()

    def _random(self,a,b,n):
        return a+(b-a)*np.random.random(n)

    def _gaussian(self,x, h, mu, sig):
        return h*np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

    def _gendata(self,nspec,nchan,ndet):
        data = np.zeros((nspec,nchan,ndet),dtype=np.float32)

        npeaks = 10
        x = np.arange(nchan,dtype=np.float32)
        h = self._random(1000,2000,npeaks)
        mu = self._random(0,nchan,npeaks)
        sig = self._random(20,30,npeaks)

        detg = np.linspace(1,1.5,ndet)

        for i in range(nspec):
            for k in range(npeaks):
                data[i,:,0] += self._gaussian(x,h[k],mu[k],sig[k])

            for j in range(ndet):
                data[i,:,j] = detg[j]#*data[i,:,0]

        return data + np.random.poisson(data)

    def _genstats(self,nspec,ndet):
        stats = np.zeros((nspec,6,ndet),dtype=np.int32)
        for i in range(ndet):
            stats[...,i] = i 
        return stats

    def _genline(self,radix,mapnum,linenum,ndet,nspec,nchan):
        line = xiaedfnd.xialine(self.dir.path,radix,mapnum,linenum,overwrite=True)

        data = self._gendata(nspec,nchan,ndet)
        stats = self._genstats(nspec,ndet)

        xialabels = ["{:02d}".format(i) for i in range(ndet)]

        line.save(data,xialabels,stats=stats)
    
        self.dir.compare(tuple(["test_xia{:02}_{:04}_0000_{:04}.edf".format(det,mapnum,linenum) for det in range(ndet)]+\
                               ["test_xiast_{:04}_0000_{:04}.edf".format(mapnum,linenum)]),path=self.dir.path)


        np.testing.assert_array_equal(data,line.data)
        np.testing.assert_array_equal(stats,line.stats)

        if ndet>1:
            line.skiponread([0])
            np.testing.assert_array_equal(data[...,1:],line.data)
            np.testing.assert_array_equal(stats[...,1:],line.stats)

        if ndet>2:
            line.skiponread([0,1])
            np.testing.assert_array_equal(data[...,2:],line.data)
            np.testing.assert_array_equal(stats[...,2:],line.stats)

        #import matplotlib.pyplot as plt
        #plt.plot(data[0,:,0])
        #plt.show()

    def test_line(self):
        for ndet in range(1,8,4):
            for nchan in [1,1024,2047]:
                for nspec in [1,10,100]:
                    self._genline("test",2,10,ndet,nspec,nchan)

def test_suite_all():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_xiaedfnd("test_line"))
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite_all()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
