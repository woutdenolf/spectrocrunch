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

from ..classfactory import FactoryMeta
from future.utils import with_metaclass

class test_classfactory(unittest.TestCase):

    def test_inheritance(self):
        class root(object):
            def __init__(self,arg0,kw0=None):
                self.arg0 = arg0
                self.kw0 = kw0
                super(root,self).__init__()
    
        class class1(with_metaclass(FactoryMeta,root)):
            def __init__(self,arg0,arg1,kw1=None,**kwargs):
                self.arg1 = arg1
                self.kw1 = kw1
                super(class1, self).__init__(arg0,**kwargs)
        factory1 = class1.factory
        
        class class2(with_metaclass(FactoryMeta,root)):
            def __init__(self,arg0,arg2,kw2=None,**kwargs):
                self.arg2 = arg2
                self.kw2 = kw2
                super(class2, self).__init__(arg0,**kwargs)
        factory2 = class2.factory
        
        class class1a(class1):
            aliases = ["class 1a"]
            def __init__(self,arg0,arg1,arg1a,kw1a=None,**kwargs):
                self.arg1a = arg1a
                self.kw1a = kw1a
                super(class1a, self).__init__(arg0,arg1,**kwargs)
                
        class class1b(class1):
            aliases = ["class 1b"]
            def __init__(self,arg0,arg1,arg1b,kw1b=None,**kwargs):
                self.arg1b = arg1b
                self.kw1b = kw1b
                super(class1b, self).__init__(arg0,arg1,**kwargs)

        class class12a(class1,class2):
            aliases = ["class 12a"]
            def __init__(self,arg0,arg1,arg2,arg12a,kw12a=None,**kwargs):
                self.arg12a = arg12a
                self.kw12a = kw12a
                super(class12a, self).__init__(arg0=arg0,arg1=arg1,arg2=arg2,**kwargs)
        
        self.assertEqual(class1.clsregistry.keys(),['class1', 'class1a', 'class1b', 'class12a'])
        self.assertEqual(class2.clsregistry.keys(),['class2', 'class12a'])
        
        o = factory1("class1","0","1",kw0="k0",kw1="k1")
        self.assertEqual((o.arg0,o.arg1,o.kw0,o.kw1),("0","1","k0","k1"))
        self.assertEqual(class1,o.__class__)
        
        o = factory1("class1a","0","1","a",kw0="k0",kw1="k1",kw1a="k1a")
        self.assertEqual((o.arg0,o.arg1,o.arg1a,o.kw0,o.kw1,o.kw1a),("0","1","a","k0","k1","k1a"))
        self.assertEqual(class1a,o.__class__)
        
        o = factory1("class 1a","0","1","a",kw0="k0",kw1="k1",kw1a="k1a")
        self.assertEqual((o.arg0,o.arg1,o.arg1a,o.kw0,o.kw1,o.kw1a),("0","1","a","k0","k1","k1a"))
        self.assertEqual(class1a,o.__class__)
        
        o = factory1("class1b","0","1","b",kw0="k0",kw1="k1",kw1b="k1b")
        self.assertEqual((o.arg0,o.arg1,o.arg1b,o.kw0,o.kw1,o.kw1b),("0","1","b","k0","k1","k1b"))
        self.assertEqual(class1b,o.__class__)
        
        o = factory1("class 1b","0","1","b",kw0="k0",kw1="k1",kw1b="k1b")
        self.assertEqual((o.arg0,o.arg1,o.arg1b,o.kw0,o.kw1,o.kw1b),("0","1","b","k0","k1","k1b"))
        self.assertEqual(class1b,o.__class__)
        
        o = factory1("class12a","0","1","2","a",kw0="k0",kw1="k1",kw2="k2",kw12a="k12a")
        self.assertEqual((o.arg0,o.arg1,o.arg2,o.arg12a,o.kw0,o.kw1,o.kw2,o.kw12a),("0","1","2","a","k0","k1","k2","k12a"))
        self.assertEqual(class12a,o.__class__)
        
        o = factory1("class 12a","0","1","2","a",kw0="k0",kw1="k1",kw2="k2",kw12a="k12a")
        self.assertEqual((o.arg0,o.arg1,o.arg2,o.arg12a,o.kw0,o.kw1,o.kw2,o.kw12a),("0","1","2","a","k0","k1","k2","k12a"))
        self.assertEqual(class12a,o.__class__)
        
        o = factory2("class2","0","2",kw0="k0",kw2="k2")
        self.assertEqual((o.arg0,o.arg2,o.kw0,o.kw2),("0","2","k0","k2"))
        self.assertEqual(class2,o.__class__)
        
        o = factory2("class12a","0","1","2","a",kw0="k0",kw1="k1",kw2="k2",kw12a="k12a")
        self.assertEqual((o.arg0,o.arg1,o.arg2,o.arg12a,o.kw0,o.kw1,o.kw2,o.kw12a),("0","1","2","a","k0","k1","k2","k12a"))
        self.assertEqual(class12a,o.__class__)
        
        o = factory2("class 12a","0","1","2","a",kw0="k0",kw1="k1",kw2="k2",kw12a="k12a")
        self.assertEqual((o.arg0,o.arg1,o.arg2,o.arg12a,o.kw0,o.kw1,o.kw2,o.kw12a),("0","1","2","a","k0","k1","k2","k12a"))
        self.assertEqual(class12a,o.__class__)
        
def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()

    testSuite.addTest(test_classfactory("test_inheritance"))
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
