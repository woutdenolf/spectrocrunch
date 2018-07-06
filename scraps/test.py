# -*- coding: utf-8 -*-
#
#   Copyright (C) 2015 European Synchrotron Radiation Facility, Grenoble, France
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

import os, sys
sys.path.insert(1,os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

#from spectrocrunch.align.tests.test_teststack import test_suite_all
#import unittest

import re

if __name__ == '__main__':
#    mysuite = test_suite_all()
#    runner = unittest.TextTestRunner()
#    if not runner.run(mysuite).wasSuccessful():
#        sys.exit(1)

    title = "scan 0  zapimage  sampy 23.962 74.962 255 100 sampz 28.252 70.452 211 0  date : Sun Sep 13 02:24:13 2015;"

    fnumber = "(?:[+-]?[0-9]*\.?[0-9]+)"
    inumber = "\d+"
    blanks = "\s+"
    motor = "[a-zA-Z]+"
    expr = "zapimage" + blanks +\
           "("+ motor +")" + blanks +\
           "("+ fnumber +")" + blanks +\
           "("+ fnumber +")" + blanks +\
           "("+ inumber +")" + blanks +\
           "("+ inumber +")" + blanks +\
           "("+ motor +")" + blanks +\
           "("+ fnumber +")" + blanks +\
           "("+ fnumber +")" + blanks +\
           "("+ inumber +")"


    print(re.findall(expr,title))
