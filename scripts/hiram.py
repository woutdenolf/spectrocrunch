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

from spectrocrunch.process.id21_fluoxas import process

if __name__ == '__main__':

    # Data
    scanname = "fXANES5"
    scannumbers = range(1,173) # from a to b-1 !!!!!
    sourcepath =  "/data/id21/inhouse/15apr/Hiram/fXANES5"

    # Fitting
    cfgfile = "/data/id21/inhouse/15apr/Hiram/results/fXANES5.cfg" # or None

    # Results
    destpath = "/data/id21/inhouse/15apr/Hiram/results/fXANES5"

    skippreprocessing = True

    # Alignment
    alignmethod = "fft" #None, fft, sift, elastix
    alignreference = "Ca-K"

    # Do processing
    process(sourcepath,destpath,scanname,scannumbers,cfgfile,alignmethod,alignreference,skippre=skippreprocessing)

