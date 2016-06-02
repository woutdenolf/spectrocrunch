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

from spectrocrunch.visualization.id21_scanoverlap import plot
import numpy as np

if __name__ == '__main__':

    specfilename = "/data/id21/store/backup_visitor/2016/hg79/id21/spec/16051103.dat"
    specnumbers = range(5,110,4)#+range(112,152,3)
    hdf5filename = "/data/id21/inhouse/xrfalign/hg79/virginsphere/monnight_maps.align.h5"
    offsamz = -34
    offsamy = -40

    #specfilename = "/data/id21/store/backup_visitor/2016/hg79/id21/spec/16050901.dat"
    #specnumbers = range(25,118,4)+range(123,169,3)+range(172,191,3)+range(194,207,4)
    #hdf5filename = "/data/id21/inhouse/xrfalign/hg79/blacksphere/RSM_BlackSphere.align.h5"
    #offsamz = -55
    #offsamy = -50

    grps = {}
    i = 0
    grps[i] = {"path":"/detector0/Pb-M","ind":0,"lo":0.05,"hi":0.20}
    i += 1
    grps[i] = {"path":"/detector0/Pb-M","ind":1,"lo":0.05,"hi":0.20}
    i += 1
    grps[i] = {"path":"/detector0/Cl-K","ind":2,"lo":0.05,"hi":0.20}
    i += 1

    grps = {}
    i = 0
    grps[i] = {"path":"/detector0/Pb-M","ind":3,"lo":0.05,"hi":0.00}
    i += 1
    grps[i] = {"path":"/detector0/Cl-K","ind":3,"lo":0.05,"hi":0.50}
    i += 1
    grps[i] = {"path":"/detector0/S-K","ind":3,"lo":0.05,"hi":0.50}
    i += 1
    specnumbers = []

    plot(hdf5filename,grps,specfilename,specnumbers,offsamy,offsamz)

