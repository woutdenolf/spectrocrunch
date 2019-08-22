# -*- coding: utf-8 -*-

import os, sys
sys.path.insert(1,os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

#from spectrocrunch.align.tests.test_teststack import test_suite
#import unittest

import re

if __name__ == '__main__':
#    mysuite = test_suite()
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
