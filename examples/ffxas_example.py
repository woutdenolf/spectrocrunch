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

# Don't use the installed version
import os, sys
sys.path.insert(1,os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import logging
logging.getLogger("spectrocrunch").setLevel(logging.INFO)

from spectrocrunch.process.id21_ffxas import process

if __name__ == '__main__':

    path = os.path.dirname(os.path.abspath(__file__))

    # Raw data
    sourcepath = "/data/id21/store/backup_visitor/2015/hg32/id21/ff/data/default/default4/ScreamSample3ab/"
    radix = "ScreamSample3ab"
    rebin = (1,1)
    skippreprocessing = True

    # Result
    destpath = os.path.join(path,"testresults",radix)

    # Normalization
    skipnormalization = True

    # Alignment
    alignmethod = "fft"
    refimageindex = 10
    crop = True
    roi = ((100,260),(160,350))
    plot = True

    # Process
    process(sourcepath,destpath,radix,rebin,alignmethod,skippre=skippreprocessing,skipnormalization=skipnormalization,\
            refimageindex=refimageindex,crop=crop,roi=roi,plot=plot)


