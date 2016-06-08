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

# Don't use the installed version
import os, sys
sys.path.insert(1,os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import logging
logging.getLogger("spectrocrunch").setLevel(logging.INFO)

from spectrocrunch.process.id21_fluoxas import process
from spectrocrunch.process.id21_fluoxas import defaultstack

if __name__ == '__main__':
    #defaultstack('/mntdirect/_data_id21_inhouse/wout/dev/SpectroCrunch/examples/testresults/fXANES5/fXANES5.norm.h5',['/detector0/Ca-K'],'Ca-K')
    #exit()

    path = os.path.dirname(os.path.abspath(__file__))

    example = "example1"

    cfgfile = None
    skippreprocessing = False
    skipnormalization = False
    dtcor = False
    default = None
    crop = False
    roi = None
    plot = True

    if example=="example1":
        scanname = ["fXANES5"]
        scannumbers = [range(1,173)]
        sourcepath = [os.path.join(path,"testdata","xrfxanes","id21",scanname[0])]
        destpath = os.path.join(path,"testresults",scanname[0])
        cfgfile = os.path.join(path,"testdata","xrfxanes","id21",scanname[0]+".cfg")

        skippreprocessing = True
        skipnormalization = False

        alignmethod = "max" #None, fft, sift, elastix
        alignreference = "Ca-K"
        refimageindex = 0 # None for pair-wise alignment
        crop = True
        roi = ((10,60),(10,60))

        default = "Ca-K"

    elif example=="example2":
        scanname = ["saphir"]
        scannumbers = [range(1,143)]
        sourcepath = [os.path.join(path,"testdata","xrfxanes","id21",scanname[0])]
        destpath = os.path.join(path,"testresults",scanname[0])
        cfgfile = os.path.join(path,"testdata","xrfxanes","id21",scanname[0]+".cfg")
        scanname[0] += "1_xrf_ff"

        skippreprocessing = True
        skipnormalization = False

        alignmethod = None #None, fft, sift, elastix
        alignreference = "xmap_x3c"
        refimageindex = 0 # None for pair-wise alignment

    elif example=="example3":
        scanname = ["roi4_p1","roi4_p1"]
        scannumbers = [range(1,54),range(1,63-10)]
        sourcepath = ["/data/visitor/in946/id21/fXAS_sample1/roi4_p1_fluoXAS_7",
                      "/data/visitor/in946/id21/fXAS_sample1/roi4_p1_fluoXAS_8"]
        destpath = os.path.join(path,"testresults","roi4_p1")
        cfgfile = None

        skippreprocessing = False
        skipnormalization = False

        alignmethod = "elastix" #None, fft, sift, elastix
        alignreference = "xmap_x3c"
        refimageindex = 6 # None for pair-wise alignment

    elif example=="example4":
        scanname = ["fXANESpl"]
        scannumbers = [range(1,5)+range(50,53)]#[range(1,84)+range(85,130)] # from a to b-1 !!!!!
        sourcepath = ["/data/id21/inhouse/15fev/Hiram/fXANESplasma"]

        # Fitting
        cfgfile = None#"/data/id21/inhouse/15apr/Hiram/results/fXANES5.cfg" # or None

        # Results
        destpath = os.path.join(path,"testresults",scanname[0])

        skippreprocessing = False
        skipnormalization = False

        # Alignment
        alignmethod = "max" #None, "fft", "sift", "elastix"
        alignreference = "xmap_x3c"
        refimageindex = 1 # None or a number

    else:
        sys.exit()

    process(sourcepath,destpath,scanname,scannumbers,cfgfile,alignmethod,alignreference,default=default, \
            refimageindex=refimageindex,skippre=skippreprocessing,skipnormalization=skipnormalization,dtcor=dtcor,crop=crop,roi=roi,plot=plot)

