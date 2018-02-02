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
import json

from .. import units
from .. import persistence

class TestPersistentClass(object):

    def __init__(self,x,y,z=None):
        self.x = x
        self.y = y
        self.z = z

    def __eq__(self,other):
        return self.x==other.x and self.y==other.y and self.z==other.z

    def serialize(self):
        x = self.x
        y = self.y
        z = self.z
        if hasattr(x,"serialize"):
            x = x.serialize()
        if hasattr(y,"serialize"):
            y = y.serialize()
        if hasattr(z,"serialize"):
            z = z.serialize()

        return persistence.SerializedGenerator(module=__name__,\
                generator="TestPersistentClass",\
                args=(x,y),\
                kwargs={"z":z})


class test_persistence(unittest.TestCase):

    def setUp(self):
        self.dir = TempDirectory()

    def tearDown(self):
        self.dir.cleanup()
        
    def test_units(self):
        a = units.Quantity(5,"keV")
        b = 10.
        c = TestPersistentClass(10,20,z=30)
        d = TestPersistentClass([1,2,3],\
                TestPersistentClass(-10,-20,z=-30),\
                z=TestPersistentClass(-10,TestPersistentClass(-10,-20,z=-30),z=TestPersistentClass(-10,-20,z=-30)))
    
        jsonfile = os.path.join(self.dir.path,"test.json")
        data = {"a":a,"b":b,"c":c,"d":d}
        
        serialized = {k:persistence.serialize(v) for k,v in data.items()}
        
        with open(jsonfile,'w') as f:
            json.dump(serialized,f,indent=2)
        
        deserialized = persistence.deserialize(jsonfile)
        
        self.assertEqual(data,deserialized)
        
def test_suite_all():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()

    testSuite.addTest(test_persistence("test_units"))
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite_all()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
