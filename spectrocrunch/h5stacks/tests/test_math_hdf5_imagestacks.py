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
import os
import numpy as np
import h5py

from .. import math_hdf5_imagestacks as math
import ...io.nexus as nexus

class test_math_hdf5_imagestacks(unittest.TestCase):
    def setUp(self):
        self.dir = TempDirectory()

    def tearDown(self):
        self.dir.cleanup()

    def test_expression(self):
        nh = 1
        nv = 1
        nd = 1

        a = np.random.rand(nv,nh,nd)
        b = np.random.rand(nv,nh,nd)

        x = np.random.rand(nv,nh,nd)
        y = np.random.rand(nv,nh,nd)
        u = np.random.rand(nv,nh,nd)
        v = np.random.rand(nv,nh,nd)
        w = np.random.rand(nv,nh,nd)

        stacka = "/stacks/a"
        stackb = "/stacks/b"
        stackx = "/stacks/x"
        stacky = "/stacks/y"
        stacku = "/stacks/u"
        stackv = "/stacks/v"
        stackw = "/stacks/w"

        innames = {"a":stacka,"b":stackb}
        args = {"a":a,"b":b,"x":x,"y":y,"u":u,"v":v,"w":w}
        stackaxes = [{"name":"v","data":np.arange(nv)},\
                     {"name":"h","data":np.arange(nh)},\
                     {"name":"d","data":np.arange(nd)}]

        expression = "{}*({x}+2*{y})/({}/3.+{u}+2*{v}+2*{w})"
        sol = [p*(x+2*y)/(p/3.+u+2*v+2*w) for p in [a,b]]
        self.evalexpression(expression,stackaxes,args,innames,sol)

    def evalexpression(self,expression,stackaxes,args,innames,sol):
        fin = os.path.join(self.dir.path,"in.h5")
        fout = os.path.join(self.dir.path,"out.h5")

        with h5py.File(fin,'w') as f:
            axes = nexus.createaxes(f,stackaxes)
            grp = nexus.newNXentry(f,"stacks")
            for k,v in args.items():
                nxdatagrp = nexus.newNXdata(grp,k,"")
                dset = nexus.createNXdataSignal(nxdatagrp,data = v)
                nexus.linkaxes(f,axes,[nxdatagrp])
                expression = math.replaceexpression(expression,[k],["/stacks/{}".format(k)])

        stacks_out, axes_out = math.calc_hdf5_imagestacks(fin,fout,axes,expression,innames.values(),innames.values(),overwrite=True)

        with h5py.File(fout,'r') as f:
            for name,expect in zip(innames,sol):
                grp = f[innames[name]]
                calc = grp[grp.attrs["signal"]]
                np.testing.assert_almost_equal(calc,expect)
                

def test_suite_all():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_math_hdf5_imagestacks("test_expression"))
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite_all()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
