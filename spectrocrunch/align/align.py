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
    """Allows for the alignment of several stacks based on one stack.
       "Alignment" is the process of determining a transformation between
       two images that should represent the same thing but transformed/deformed.
    """

    def __init__(self,source,sourcelist,dest,destlist,extension,stackdim=None,overwrite=False,cval=np.nan,plot=False):

        # Data IO
        self.source = alignSource(source,sourcelist,stackdim=stackdim)
        self.dest = alignDest(dest,destlist,extension,stackdim=stackdim,overwrite=overwrite)
    
        # Missing data
        self.cval = cval #TODO: use nan's as default but pymca can't handle them

        # Transformation settings (set before actual transformation)
        self.dtype = np.float32
        self.alignonraw = True
        self.usekernel = False # Doesn't work well for Elastix!
        self.pre_align = {"ROI":None}
        self.pre_transform = {"extended_dimensions":True}
        self.padding = ((0,0),(0,0))

        # Affine transformation (change of frame matrices, not change of coordinates!)
        self.offsets = np.zeros((self.source.nimages,2),dtype=self.dtype)
        self.linears = np.zeros((self.source.nimages,2,2),dtype=self.dtype)
        self.origin = np.zeros(2,dtype = self.dtype)

        # Others
        self.idlinear = np.identity(2,dtype = self.dtype)
        self.idoffset = np.zeros(2,dtype = self.dtype)

        self.doplot = plot

    def plot(self,img,index,title):
        """Visualize alignment in progress
        """
        if not self.doplot:
            return
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
        """Apply image padding.
        """
        return np.pad(img,self.padding,'constant',constant_values=self.cval)

    def roi(self,img,roi):
        """Extract ROI
        """
        n1,n2 = img.shape
        a1 = roi[0][0]
        b1 = roi[0][1]
        a2 = roi[1][0]
        b2 = roi[1][1]
        if a1 < 0:
            a1 += n1
        if b1 < 0:
            b1 += n1
        if a2 < 0:
            a2 += n2
        if b2 < 0:
            b2 += n2
        if b1 <= a1 or b2 <= a2 or \
           a1 < 0 or a2 < 0 or b1 < 0 or b2 < 0 or \
           a1 >= n1 or a2 >= n2 or b1 >= n1 or b2 >= n2:
           raise ValueError("ROI is invalid.")
        return img[a1:b1+1,a2:b2+1]

    def writeimg(self,img,datasetindex,imageindex):
        """Save 1 image in 1 stack.
        """
        self.dest.writeimg(img,datasetindex,imageindex)

    def copyimg(self,datasetindex,imageindex):
        """Copy 1 image in 1 stack.
        """
        img = self.readimgraw(datasetindex,imageindex)
        self.writeimg(img,datasetindex,imageindex)

    def readimgraw(self,datasetindex,imageindex):
        """Read 1 image in 1 stack.
        """
        return self.source.readimgas(datasetindex,imageindex,self.dtype)

    def readimgrawprep(self,datasetindex,imageindex):
        """Get raw image, preprocessed for alignment
        """
        img = self.readimgraw(datasetindex,imageindex)
        img = self.dopre_align(img)
        return img

    def nopre_align(self):
        """
        Returns:
            bool: True when alignment is done on the raw image
        """
        return self.pre_align["ROI"] is None

    def dopre_align(self,img):
        if self.pre_align["ROI"] is not None:
            img = self.roi(img,self.pre_align["ROI"])
        return img

    def nopre_transform(self):
        """
        Returns:
            bool: True when transformation is done on the raw image
        """
        return self.pre_transform["extended_dimensions"]==False

    def dopre_transform(self,img):
        """Manual transformation before the real transformation (not used in alignment)
        """
        if self.pre_transform["extended_dimensions"]:
            img = self.pad(img)
        return img

    def execute_alignkernel(self,img):
        raise NotImplementedError()

    def execute_transformkernel(self,img):
        raise NotImplementedError()

    def execute_transform_nokernel(self,img,offset,linear):
        """Apply a transformation to an image
        """
        # affine_transform: assumes "change of frame" matrices with coordinates (y,x)
        # shift: assumes "change of coordinates" matrices with coordinates (y,x)
        if np.array_equal(linear,self.idlinear):
            return scipy.ndimage.interpolation.shift(img,-offset,cval = self.cval,order=1,mode="constant")
        else:
            raise NotImplementedError()
            # Use origin!
            return scipy.ndimage.interpolation.affine_transform(img,linear,offset=offset,cval = self.cval,order=1,mode="constant")

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
        return self.nopre_transform() and self.transformidentity(i)

    def transform(self,img,i):
        """Apply image transformation
        """
        # Return when transformation is the identity
        if self.pureidentity(i):
            return img

        # Apply initial transformation (not used in alignment)
        imgtransformed = self.dopre_transform(img)

        # Return when transformation is the identity
        if self.transformidentity(i):
            return imgtransformed

        # Apply transformation
        if self.usekernel:
            imgtransformed = self.execute_transformkernel(imgtransformed)
        else:
            imgtransformed = self.execute_transform(imgtransformed,i)

        return imgtransformed

    def calculateorigin(self):
        """Because the image has a different preprocessing for alignment and transformation,
           the origin of the image used during alignment in the image used during transformation
           must be calculated.
        """
        self.origin = np.copy(self.idoffset)

        if self.pre_align["ROI"] is not None:
            self.origin[0] += self.pre_align["ROI"][0][0]
            self.origin[1] += self.pre_align["ROI"][1][0]

        if self.pre_transform["extended_dimensions"]:
            self.origin[0] += self.padding[0][0]
            self.origin[1] += self.padding[1][0]

    def get_transformation(self):
        raise NotImplementedError()

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

    def set_transformation(self,offset,changed):
        raise NotImplementedError()

    def update_transformation(self,i,bchange):
        """Update transformation parameters
        """
        if self.usekernel:
            self.set_transformation(self.offsets[i,...],bchange)

    def settransformidentity(self,i):
        """Make this transformation the identity
        """
        self.offsets[i,...] = self.idoffset
        self.linears[i,...] = self.idlinear

    def minimaltransformation(self):
        """If all transformations are known, they can be reduced to minimize the difference with the original images
        """
        self.offsets[:,0] -= np.median(self.offsets[:,0])
        self.offsets[:,1] -= np.median(self.offsets[:,1])

    def paddingfromtransformation(self):
        """If all transformations are known, padding can be calculated so nothing goes out of the field of view
        """
        paddim1before = np.ceil(max(np.max(self.offsets[:,0]),0)).astype(np.int)
        paddim1after = np.ceil(max(-np.min(self.offsets[:,0]),0)).astype(np.int)
        paddim2before = np.ceil(max(np.max(self.offsets[:,1]),0)).astype(np.int)
        paddim2after = np.ceil(max(-np.min(self.offsets[:,1]),0)).astype(np.int)
        self.padding = ((paddim1before,paddim1after),(paddim2before,paddim2after))

    def parsetransformation_beforeapplication(self):
        """Adapt transformations before applying them
        """
        self.minimaltransformation()
        self.paddingfromtransformation()

    def getaxesaftertransformation(self,axes):
        """Image X and Y axes after transformation
        """
        if not self.pre_transform["extended_dimensions"]:
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
        """Allocate space for saving results
        """
        nimages = self.source.nimages
        imgsize = self.source.imgsize
        
        if self.pre_transform["extended_dimensions"]:
            imgsize = (imgsize[0] + self.padding[0][0] + self.padding[0][1],
                       imgsize[1] + self.padding[1][0] + self.padding[1][1])
            
        self.dest.prepare(nimages,imgsize,self.dtype)

    def set_reference(self,img,previous=False):
        raise NotImplementedError()

    def doalign(self,refdatasetindex,refimageindex=None,aligntype=alignType.full):
        """Align datasets and save the result.
           Calculate transformation (alignment):
            fixed image: raw -> align prep
            moving image: raw -> align prep -> transformed
           Apply transformation:
            aligned image: raw -> transform prep -> transformed
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
                    self.settransformidentity(i)
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
                    if j==refdatasetindex and self.nopre_align() and self.nopre_transform() and aligntype!=alignType.usetransfo:
                        img = imgaligned
                    else:
                        img = self.readimgraw(j,i)
                        img = self.transform(img,i)

                    self.writeimg(img,j,i)

    def align(self,refdatasetindex,refimageindex = None,onraw = False,extend = False,redo = False,roi = None):
        """Alignment function that needs to be called 
        
        Args:
            refdatasetindex(int): stack to be used for alignment
            refimageindex(Optional(int)): fixed reference to align on (pairwise alignment when None)
            onraw(Optional(bool)): when doing pairwise alignment, use the previous raw or aligned images to align the next one on
            extend(Optional(bool)): make sure nothing is transformed outside the field of view
            redo(Optional(bool)): apply transformations without recalculating them
        """

        pairwise = refimageindex is None
        if pairwise:
            self.alignonraw = onraw
        else:
            self.alignonraw = True

        self.pre_align["ROI"] = roi
        self.pre_transform["extended_dimensions"] = extend
        self.calculateorigin()

        if redo:
            self.doalign(refdatasetindex,aligntype=alignType.usetransfo)
        else:
            if extend:
                self.doalign(refdatasetindex,refimageindex=refimageindex,aligntype=alignType.calctransfo)
                self.parsetransformation_beforeapplication()
                self.doalign(refdatasetindex,refimageindex=refimageindex,aligntype=alignType.usetransfo)
            else:
                self.doalign(refdatasetindex,refimageindex=refimageindex)

