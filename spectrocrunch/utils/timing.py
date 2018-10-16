# -*- coding: utf-8 -*-
#
#   Copyright (C) 2016 European Synchrotron Radiation Facility, Grenoble, France
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

import datetime
import sys
import contextlib
import time

def taketimestamp():
    return datetime.datetime.now()

def hms(seconds):
    sec = seconds
    hours = int(sec/3600.)
    sec -= hours*3600
    min = int(sec/60.)
    sec -= min*60
    sec = int(sec)
    return (hours,min,sec)

def strseconds(seconds):
    return "{:d}h {:d}m {:d}s".format(*hms(seconds))

def printtimeelapsed(T0,logger,text="Elapsed time"):
    hours,min,sec = hms((datetime.datetime.now()-T0).seconds)
    logger.info("  %s: %dh %dm %ds" % (text,hours,min,sec))

class progress(object):
    def __init__(self,logger):
        self.logger = logger
        self.setn(1)
        self.start()

    def start(self):
        self.T0 = datetime.datetime.now()

    def setn(self,nmax):
        self.nmax = nmax
        self.ncurrent = 0
        self.setnfine(1)

    def setnfine(self,nmaxfine):
        self.nmaxfine = nmaxfine
        self.ncurrentfine = 0

    def ndone(self,n):
        self.ncurrent += n
        self.ncurrentfine = 0

    def ndonefine(self,n):
        self.ncurrentfine += n

    def printprogress(self):
        # Ttot = Telapsed + Tleft
        # Telapsed/Ttot = ncurrent/nmax = pdone
        # Tleft = Ttot - Telapsed = Telapsed*(1/pdone - 1) 

        nmax = self.nmax*self.nmaxfine
        ncurrent = self.ncurrent*self.nmaxfine + self.ncurrentfine

        Telapsed = (datetime.datetime.now()-self.T0).seconds
        
        hours,min,sec = hms(Telapsed)

        hours2,min2,sec2 = hms(Telapsed*(nmax/float(ncurrent)-1))

        pdone = int(ncurrent/float(nmax)*100)

        self.logger.info("  Elapsed: %dh %dm %ds (%d%%)  Left: %dh %dm %ds (%d%%)" % (hours,min,sec,pdone,hours2,min2,sec2,100-pdone))
        sys.stdout.flush()

@contextlib.contextmanager
def timeit(name=""):
    t0 = time.time()
    yield
    t1 = time.time()

    print("Execution time ({}): {}".format(name,t1-t0))
    
