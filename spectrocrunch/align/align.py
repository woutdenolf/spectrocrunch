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

import scipy.ndimage.interpolation
import numpy as np
import pylab

from .alignSource import alignSource
from .alignDest import alignDest
from .types import alignType

class align(object):

    def __init__(self,source,sourcelist,dest,destlist,extension,stackdim=None,overwrite=False,cval=0):
        # Data IO
        self.source = alignSource(source,sourcelist,stackdim=stackdim)
        self.dest = alignDest(dest,destlist,extension,stackdim=stackdim,overwrite=overwrite)
    
        # Missing data
        self.cval = cval #TODO: use nan's as default but pymca can't handle them

        # Transformation settings (set before actual transformation)
        self.dtype = np.float32
        self.alignonraw = True
        self.usekernel = False
        self.extended_dimensions_alignment = False
        self.extended_dimensions_initialtransform = True
        self.padding = ((0,0),(0,0))

        # Affine transformation (change of frame matrices, not change of coordinates!)
        self.offsets = np.zeros((self.source.nimages,2),dtype=self.dtype)
        self.linears = np.zeros((self.source.nimages,2,2),dtype=self.dtype)

        # Others
        self.idlinear = np.identity(2,dtype = self.dtype)
        self.idoffset = np.zeros(2,dtype = self.dtype)

    def plot(self,img,index,title):
        #pylab.figure(index)
        #pylab.clf()

        pylab.figure(1)
        pylab.subplot(130+index)
        img2 = img.copy()
        img2[np.isnan(img2)] = 0
        pylab.imshow(img2,origin='lower',interpolation='nearest')
        pylab.title(title)
        pylab.pause(0.01)

    def pad(self,img):
        return np.pad(img,self.padding,'constant',constant_values=self.cval)

    def readimgraw(self,datasetindex,imageindex):
        """Get raw image
        """
        return self.source.readimgas(datasetindex,imageindex,self.dtype)

    def readimgrawprep(self,datasetindex,imageindex):
        """Get raw image, preprocessed for alignment
        """
        img = self.readimgraw(datasetindex,imageindex)

        if self.extended_dimensions_alignment:
            img = self.pad(img)

        return img

    def writeimg(self,img,datasetindex,imageindex):
        self.dest.writeimg(img,datasetindex,imageindex)

    def copyimg(self,datasetindex,imageindex):
        img = self.readimgraw(datasetindex,imageindex)
        self.writeimg(img,datasetindex,imageindex)

    def alignprepidentity(self):
        """Is there any raw preprocessing?
        """
        return self.extended_dimensions_alignment == False

    def initialtransform(self,img):
        """Manual transformation before the real transformation (not used in alignment)
        """
        if self.extended_dimensions_initialtransform:
            img = self.pad(img)
        return img
    
    def initialidentity(self):
        """Is the transformation the identity
        """
        return self.extended_dimensions_initialtransform==False

    def execute_transform_nokernel(self,img,offset,linear):
        # affine_transform: assumes "change of frame" matrices with coordinates (y,x)
        # shift: assumes "change of coordinates" matrices with coordinates (y,x)
        if np.array_equal(linear,self.idlinear):
            return scipy.ndimage.interpolation.shift(img,-offset,cval = self.cval,order=1,mode="constant")
        else:
            return scipy.ndimage.interpolation.affine_transform(data,linear,offset=offset,cval = self.cval,order=1,mode="constant")

    def execute_transform(self,img,i):
        """Transform according to the transformation extracted from an the transformation kernel (see gettransformation).
        """
        return self.execute_transform_nokernel(img,self.offsets[i,...],self.linears[i,...])

    def transformidentity(self,i):
        """Is the transformation the identity
        """
        return np.array_equal(self.offsets[i,...], self.idoffset) and np.array_equal(self.linears[i,...], self.idlinear)

    def pureidentity(self,i):
        """Is the transformation the identity, including the changes applied before (like padding)
        """
        return self.initialidentity() and self.transformidentity(i)

    def transform(self,img,i):
        """Like transform but including the changes applied before (like padding)
        """
        # Return when transformation is the identity
        if self.pureidentity(i):
            return img

        # Apply initial transformation (not used in alignment)
        imgtransformed = self.initialtransform(img)

        # Return when transformation is the identity
        if self.transformidentity(i):
            return imgtransformed

        # Apply transformation
        if self.usekernel:
            imgtransformed = self.execute_transformkernel(imgtransformed)
        else:
            imgtransformed = self.execute_transform(imgtransformed,i)

        return imgtransformed

    def gettransformation(self,i,pairwise):
        """Get transformation parameters
        """
        offset = self.get_transformation()

        bchange = pairwise and self.alignonraw and i!=0

        if bchange:
            self.offsets[i,...] = self.offsets[i-1,...] + offset
        else:
            self.offsets[i,...] = offset

        self.linears[i,...] = self.idlinear

        self.update_transformation(i,bchange)

    def update_transformation(self,i,bchange):
        """Update transformation parameters
        """
        if self.usekernel:
            self.set_transformation(self.offsets[i,...],bchange)

    def setaligntransformidentity(self,i):
        """Make this transformation the identity
        """
        self.offsets[i,...] = self.idoffset
        self.linears[i,...] = self.idlinear

    def minimaltransformation(self):
        self.offsets[:,0] -= np.median(self.offsets[:,0])
        self.offsets[:,1] -= np.median(self.offsets[:,1])

    def paddingfromtransformation(self):
        paddim1before = np.ceil(max(np.max(self.offsets[:,0]),0)).astype(np.int)
        paddim1after = np.ceil(max(-np.min(self.offsets[:,0]),0)).astype(np.int)
        paddim2before = np.ceil(max(np.max(self.offsets[:,1]),0)).astype(np.int)
        paddim2after = np.ceil(max(-np.min(self.offsets[:,1]),0)).astype(np.int)
        self.padding = ((paddim1before,paddim1after),(paddim2before,paddim2after))

    def parsetransformation_beforeapplication(self):
        self.minimaltransformation()
        self.paddingfromtransformation()

    def getaxesafteralignment(self,axes):
        if not self.extended_dimensions_initialtransform:
            return

        if self.source.stackdim==2:
            ind = [0,1]
        elif self.source.stackdim==1:
            ind = [0,2]
        else:
            ind = [1,2]

        for i in range(len(ind)):
            j = ind[i]

            nleft = self.padding[i][0]
            nright = self.padding[i][1]
            naxis = len(axes[j])
            newaxis = np.empty(naxis + nleft + nright,dtype = axes[j].dtype)

            off0 = 0
            axis0 = 0
            axisn = naxis
            if nleft < 0:
                axis0 -= nleft
            else:
                off0 += nleft
            if nright < 0:
                axisn += nright
            offn = off0 + axisn - axis0

            if nleft > 0:
                delta = axes[j][1]-axes[j][0]
                newaxis[0:nleft] = (axes[j][0] - delta*nleft) + delta*np.arange(nleft)

            newaxis[off0:offn] = axes[j][axis0:axisn]

            if nright > 0:
                delta = axes[j][-1]-axes[j][-2]
                newaxis[offn:] = (axes[j][-1] + delta) + delta*np.arange(nright)

            axes[j] = newaxis

    def preparedestination(self):
        nimages = self.source.nimages
        imgsize = self.source.imgsize
        
        if self.extended_dimensions_initialtransform:
            imgsize = (imgsize[0] + self.padding[0][0] + self.padding[0][1],
                       imgsize[1] + self.padding[1][0] + self.padding[1][1])

        self.dest.prepare(nimages,imgsize,self.dtype)

    def align(self,refdatasetindex,refimageindex = None,onraw = False,extend = False,redo = False):
        """Alignment function that needs to be called 
        """
        pairwise = refimageindex is None
        if pairwise:
            self.alignonraw = onraw
        else:
            self.alignonraw = True

        self.extended_dimensions_alignment = extend and not self.alignonraw
        self.extended_dimensions_initialtransform = extend

        if redo:
            self.doalign(refdatasetindex,aligntype=alignType.usetransfo)
        else:
            if extend:
                self.doalign(refdatasetindex,refimageindex=refimageindex,aligntype=alignType.calctransfo)
                self.parsetransformation_beforeapplication()
                self.doalign(refdatasetindex,refimageindex=refimageindex,aligntype=alignType.usetransfo)
            else:
                self.doalign(refdatasetindex,refimageindex=refimageindex)

    def doalign(self,refdatasetindex,refimageindex=None,aligntype=alignType.full):
        """Align datasets and save the result.
           Calculate transformation:
            fixed image: raw -> align prep
            moving image: raw -> align prep -> aligntransform
           Apply transformation:
            aligned image: raw -> initialtransform -> aligntransform
        """
        pairwise = refimageindex is None

        # Prepare destination
        if aligntype!=alignType.calctransfo:
            self.preparedestination()

        # First reference image
        if aligntype!=alignType.usetransfo:
            if pairwise:
                # Pair-wise alignment: first image is the first reference
                imgref = self.readimgrawprep(refdatasetindex,0)
                iref = 0
                self.plot(imgref,1,"Image %d (pair-wise)"%iref)
            else:
                # Fixed-reference alignment
                rawprep = self.readimgrawprep(refdatasetindex,refimageindex)
                iref = refimageindex
                self.plot(rawprep,1,"Reference %d (fixed)"%iref)
                self.set_reference(rawprep)

        # Loop over the images
        for i in range(self.source.nimages):
            if aligntype!=alignType.usetransfo:
                # Image i
                rawprep = self.readimgrawprep(refdatasetindex,i)

                # Update fixed image
                if pairwise:
                    self.set_reference(imgref)

                # Get align transformation
                if i == iref:
                    self.setaligntransformidentity(i)
                    imgaligned = rawprep
                else:
                    # Align image i to reference
                    self.plot(rawprep,3,"To align %d"%i)
                    imgaligned = self.execute_alignkernel(rawprep)
                    self.plot(imgaligned,2,"Aligned %d"%i)
                    self.gettransformation(i,pairwise)

                # Reference for the next image
                if pairwise:
                    if self.alignonraw:
                        imgref = rawprep
                    else:
                        imgref = imgaligned
                    iref = i

                # Only need transformation, not the aligned result
                if aligntype==alignType.calctransfo:
                    continue

            # Save the transformed image i of all datasets
            if self.pureidentity(i):
                for j in range(self.source.nsets):
                    self.copyimg(j,i)
            else:
                if aligntype==alignType.usetransfo:
                    self.update_transformation(i,True)

                for j in range(self.source.nsets):
                    if j==refdatasetindex and self.alignprepidentity() and self.initialidentity() and aligntype!=alignType.usetransfo:
                        img = imgaligned
                    else:
                        img = self.readimgraw(j,i)
                        img = self.transform(img,i)

                    self.writeimg(img,j,i)

