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

import matplotlib.pyplot as plt
import numpy as np

from .. import scene
from ... import ureg

class test_scene(unittest.TestCase):

    def test_images(self):
        n0,n1 = 5,10
        img = np.arange(n0*n1).reshape(n0,n1)
        
        unit0 = ureg.mm
        unit1 = ureg.micrometer
        
        s1 = scene.Scene(unit0=unit0,unit1=unit1)
        
        s2 = scene.Scene(unit0=unit0,unit1=unit1)
        s2.transpose(True)
        #s2.flipx(increasing=True)
        s2.axlabels = ["dim0","dim1"]
        s2.cmap = plt.get_cmap('gray')
        
        o1 = scene.Image(img,lim0=s1.q0([8,8+n0-1]),lim1=s1.q1([10+n1-1,10]))
        s1.register(o1)
        s2.register(o1)
        
        p0 = sorted(o1.datarange(0,border=False))
        p1 = sorted(o1.datarange(1,border=False))
        o = scene.Polyline([p0[0],p0[1],p0[1],p0[0]],[p1[0],p1[0],p1[1],p1[1]])
        s1.register(o)
        s2.register(o)
        o.set_setting("scatter",True)

        o2 = scene.Image(img,lim0=s1.q0([-2,-2+n0-1]),lim1=s1.q1([-1,-1+n1-1]))
        s1.register(o2)
        s2.register(o2)
        o.set_setting("scatter",True)
        
        p0 = sorted(o2.datarange(0,border=False))
        p1 = sorted(o2.datarange(1,border=False))
        o = scene.Text([p0[0],p0[1],p0[1],p0[0]],[p1[0],p1[0],p1[1],p1[1]],labels=[1,2,3,4])
        s1.register(o)
        s2.register(o)
        
        f, ax = plt.subplots()
        s1.setaxes(ax)
        f, ax = plt.subplots()
        s2.setaxes(ax)
        
        # Update scene 1
        s1.update()
        
        # Shift image, axes scaling and update scene 2
        o1.lim[0] = s1.q0([9,9+n0-1])
        s2.setdatarange(0,s1.q0([0,1]))
        s2.setdatarange(1,s1.q1([0,1]))
        s2.update()
        #plt.pause(0.01)
        
        # Update scene 1
        s1.update()
        
        # Reset axes of scene 1
        f, ax = plt.subplots()
        s1.setaxes(ax)
        
        # Shift image, axes offset, different normalization and update scene 1
        o1.lim[0] = s1.q0([9,9+n0-1])
        
        s1.set_setting("cnorm",scene.ColorNorm("PowerNorm",0.1))
        s1.update()
        
        #plt.pause(0.01)
        #plt.show()
        
def test_suite_all():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_scene("test_images"))
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite_all()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
