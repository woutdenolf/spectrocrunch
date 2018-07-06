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

from spectrocrunch.align.alignElastix import alignElastix
from spectrocrunch.align.alignFFT import alignFFT
from spectrocrunch.align.alignSift import alignSift
from spectrocrunch.align.alignSimple import alignMax
from spectrocrunch.align.alignSimple import alignMin
from spectrocrunch.align.types import transformationType

from spectrocrunch.align.tests.teststack import teststack

import os
import glob
import numpy as np
import h5py
import fabio
import matplotlib.pyplot as plt

def showstack(stackData,cofs):
    import PyQt4.Qt as qt
    from PyMca5.PyMca import StackBrowser

    app = qt.QApplication([])

    w = StackBrowser.StackBrowser()
    w.setStackDataObject(stackData, index=0)

    w.graphWidget.graph.setGraphTitle("hacked!")

    nimg,ny,nx = stackData.shape
    xy = np.empty((3,4))
    
    for i in range(-1):
        xy[0,:] = [0,nx-1,0,nx-1] #x
        xy[1,:] = [0,0,ny-1,ny-1] #y
        xy[2,:] = [1,1,1,1]
        M = cofs[i,...]
        if len(M)==2:
            xy = xy - np.append(M,0)[:,np.newaxis]
        else:
            xy = np.dot(np.linalg.inv(M),xy)
        ind = [0,1,3,2,0]
        print("Polygon coordinates")
        print(xy)
        x = xy[0,ind]/xy[2,ind]
        y = xy[1,ind]/xy[2,ind]
        w.graphWidget.graph.addItem(x,y, legend='border{}'.format(i))

    w.show()
    
    app.exec_()

def testelastix(img0,img1):
    import SimpleITK as sitk
    
    fig,axes = plt.subplots(1,3)

    elastix = sitk.SimpleElastix()
    elastix.SetFixedImage(sitk.GetImageFromArray(img0))
    elastix.SetMovingImage(sitk.GetImageFromArray(img1))
    p = sitk.GetDefaultParameterMap("rigid")
    elastix.SetParameterMap(p)
    elastix.Execute()
    img2 = sitk.GetArrayFromImage(elastix.GetResultImage())

    axes[0].imshow(img0,origin='lower',interpolation='nearest',cmap='jet')
    axes[0].set_title('fixed')
    axes[1].imshow(img2,origin='lower',interpolation='nearest',cmap='jet')
    axes[1].set_title('aligned')
    axes[2].imshow(img1,origin='lower',interpolation='nearest',cmap='jet')
    axes[2].set_title('moving')
    
    plt.show()
    exit()

def alignexample(t):

    # Source data (several image stacks)
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)),"testdata")

    transfotype = transformationType.translation
    onraw = False
    stackdim = 2

    roi = None
    if t=="xrfxanes":
        source = os.path.join(path,t,"h5","rape5_XANESfull.h5")
        stackdim = 2
        sourcelist = ["DSAg","DSCu","DSIt","DSZn"]
        nstack = len(sourcelist)
        refdatasetindex = 1
        refimageindex = 0
        alignclass = alignElastix
    elif t =="fluotomo":
        source = ""
        dataset1 = sorted(glob.glob(os.path.join(path,t,"*.tif")))[0:3]
        sourcelist = [dataset1]
        nstack = 1
        refdatasetindex = 0
        refimageindex = None
        onraw = True
        alignclass = alignElastix
    elif t == "testdata":
        transfotype = transformationType.translation
        source,relcof,stackdim = teststack(transfotype,nimages=5)
        #source[0]=source[0][...,0:2]
        #source[0][...,1] = source[0][...,0]
        sourcelist = None
        nstack = len(source)
        refdatasetindex = 0
        refimageindex = 0
        alignclass = alignSift
        roi = None
        onraw = True
        roi = ((5,-5),(2,-2))
        #roi = ((0,20),(60,79))
        #alignclass = alignMin
        #roi = ((10,30),(30,50))
    else:
        return

    #testelastix(source[0][...,0],source[0][...,1])

    # Destination
    outputstack = [np.zeros(1,dtype=np.float32)]*nstack

    # Align
    o = alignclass(source,sourcelist,outputstack,None,None,stackdim=stackdim,overwrite=True,plot=True,transfotype=transfotype)
    o.align(refdatasetindex,refimageindex = refimageindex,onraw = onraw,pad = True,crop = False,roi = roi)

    if t == "testdata":
        print("Know pairwise cof:")
        print(relcof)

        print("Determined absolute cofs:")
        print(o.absolute_cofs())

    # Show result
    showstack(np.rollaxis(outputstack[0],stackdim,0),o.absolute_cofs())

    # Save aligned data
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)),"testresults","fluotomo")
    if t == "fluotomo":
        with h5py.File(os.path.join(path,"fluotomo.h5")) as f:
            f["aligned"] = outputstack[0]
        with h5py.File(os.path.join(path,"fluotomo_original.h5")) as f:
            f["stack"] = np.zeros(outputstack[0].shape)
            for i in range(len(sourcelist[0])):
                f["stack"][...,i] = fabio.open(sourcelist[0][i]).data


if __name__ == '__main__':

    #alignexample("xrfxanes")
    #alignexample("fluotomo")
    alignexample("testdata")


