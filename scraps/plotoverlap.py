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
    specnumbers = range(5,110,4)
    specnumbers = [5, 9, 13, 17, 21, 25, 29, 33, 37, 41, 45, 49, 53, 57, 61, 65, 69, 73, 77, 81, 85, 89, 93, 97, 101, 105, 109]

    hdf5filename = "/data/id21/inhouse/xrfalign/hg79/virginsphere/monnight_maps.align.h5"
    offsamz = 0
    offsamy = 0

    grps = {}
    i = 0
    grps[i] = {"path":"/detector0/Pb-M","ind":3,"lo":0.05,"hi":0.00}
    i += 1
    grps[i] = {"path":"/detector0/Cl-K","ind":3,"lo":0.05,"hi":0.50}
    i += 1
    grps[i] = {"path":"/detector0/S-K","ind":3,"lo":0.05,"hi":0.50}
    i += 1

    plot(hdf5filename,grps,specfilename,specnumbers,offsamy,offsamz,transpose=False,flipvert=True,fliphor=False,defaultorigin=False,showlabels=False,color='#ffffff',printpos=False)

