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
from .types import transformationType

class align(object):
    """Allows for the alignment of several stacks based on one stack.
       "Alignment" is the process of determining a transformation between
       two images that should represent the same thing but transformed/deformed.
    """

    def __init__(self,source,sourcelist,dest,destlist,extension,stackdim=None,overwrite=False,cval=np.nan,plot=False,transfotype=transformationType.translation):

        # Data IO
        self.source = alignSource(source,sourcelist,stackdim=stackdim)
        self.dest = alignDest(dest,destlist,extension,stackdim=stackdim,overwrite=overwrite)
    
        # Missing data
        self.cval = cval #TODO: use nan's as default but pymca can't handle them

        # Transformation settings (set before actual transformation)
        self.dtype = np.float32
        self.alignonraw = True
        self.usekernel = False # Doesn't work well for Elastix!
        self.pre_align = {"roi":None}
        self.pre_transform = {"pad":False}
        self.post_transform = {"crop":False}
        self.extend = ((0,0),(0,0)) # negative: crop, positive: pad

        # Transformation (change of frame matrices, not change of coordinates!)
        self.transfotype = transfotype
        self.cofs = np.zeros((self.source.nimages,3,3),dtype=self.dtype)
        self.origin = np.zeros(2,dtype = self.dtype)

        # Others
        self.idcof = np.identity(3,dtype = self.dtype)
        self.idlinear = np.identity(2,dtype = self.dtype)
        self.idoffset = np.zeros(2,dtype = self.dtype)
        self.idproj = np.zeros(2,dtype = self.dtype)
        self.idorigin = np.zeros(2,dtype = self.dtype)

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
        """Apply padding
        """
        pad = ((max(self.extend[0][0],0),max(self.extend[0][1],0)),\
               (max(self.extend[1][0],0),max(self.extend[1][1],0)))
        if np.count_nonzero(pad)!=0:
            return np.pad(img,pad,'constant',constant_values=self.cval)
        else:
            return img

    def crop(self,img):
        """Apply cropping
        """
        dim1,dim2 = img.shape
        crop = ((max(-self.extend[0][0],0),dim1-max(-self.extend[0][1],0)),\
                (max(-self.extend[1][0],0),dim2-max(-self.extend[1][1],0)))
        if crop[0][0]!=0 or crop[1][0]!=0 or crop[0][1]!=dim1 or crop[1][1]!=dim2:
            return img[crop[0][0]:crop[0][1],crop[1][0]:crop[1][1]]
        else:
            return img

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
        return self.pre_align["roi"] is None

    def dopre_align(self,img):
        if self.pre_align["roi"] is not None:
            img = self.roi(img,self.pre_align["roi"])
        return img

    def nopre_transform(self):
        """
        Returns:
            bool: True when transformation is done on the raw image
        """
        return np.all(np.asarray(self.extend)<=0)

    def dopre_transform(self,img):
        """Manual transformation before the real transformation (not used in alignment)
        """
        if np.any(np.asarray(self.extend)>0):
            img = self.pad(img)
        return img

    def nopost_transform(self):
        """
        Returns:
            bool: True when transformation doesn't have any post processing
        """
        return np.all(np.asarray(self.extend)>=0)

    def dopost_transform(self,img):
        """Manual transformation after the real transformation (not used in alignment)
        """
        if np.any(np.asarray(self.extend)<0):
            img = self.crop(img)
        return img

    def execute_alignkernel(self,img):
        raise NotImplementedError()

    def execute_transformkernel(self,img):
        raise NotImplementedError()

    def execute_transform_nokernel(self,img,cof):
        """Apply a transformation to an image
        """

        # cof: change-of-frame matrix for coordinates (x,y)
        nolinear = np.array_equal(cof[0:2,0:2],self.idlinear)
        noproj = np.array_equal(cof[2,0:2],self.idproj)
        if nolinear and noproj:
            # shift: assumes change-of-coordinates vector for coordinates (y,x)
            return scipy.ndimage.interpolation.shift(img,-cof[1::-1,2],cval = self.cval,order=1,mode="constant")
        elif noproj:
            # affine_transform: assumes change-of-frame matrix for coordinates (y,x)
            M = self.cofwithorigin(cof)
            return scipy.ndimage.interpolation.affine_transform(img,M[0:2,0:2].T,offset=M[1::-1,2],cval = self.cval,order=1,mode="constant")
        else:
            raise NotImplementedError()

    def cofwithorigin(self,C2):
        """ C2 is calculated in a frame which is shifted with respect to the (padded/cropped) image frame.
            L3.L2.L1.X = X'
            C1.C2.C3.X' = X
            C1 = from padded/cropped frame to roi frame
            C3 = from roi frame to padded frame
        """
        
        # origin: position of the origin of the frame in which cof is calculated
        C1 = self.idcof.copy()
        C3 = self.idcof.copy()
        C1[0:2,2] = self.origin
        C3[0:2,2] = -self.origin
        return np.dot(C1,np.dot(C2,C3))

    def execute_transform(self,img,i):
        """Transform according to the transformation extracted from an the transformation kernel (see gettransformation).
        """
        return self.execute_transform_nokernel(img,self.cofs[i,...])

    def transformidentity(self,i):
        """Is the transformation the identity
        """
        return np.array_equal(self.cofs[i,...], self.idcof)

    def pureidentity(self,i):
        """Is the transformation the identity, including the changes applied before (like padding)
        """
        return self.nopre_transform() and self.nopost_transform() and self.transformidentity(i)

    def transform(self,img,i):
        """Apply image transformation
        """
        # Return when transformation is the identity
        if self.pureidentity(i):
            return img

        # Apply initial transformation (not used in alignment)
        imgtransformed = self.dopre_transform(img)

        # Apply transformation
        if not self.transformidentity(i):
            if self.usekernel:
                imgtransformed = self.execute_transformkernel(imgtransformed)
            else:
                imgtransformed = self.execute_transform(imgtransformed,i)

        # Apply final transformation (not used in alignment)
        imgtransformed = self.dopost_transform(imgtransformed)
        return imgtransformed

    def calculateorigin(self):
        """Because the image has a different preprocessing for alignment and transformation,
           the origin of the image used during alignment in the image used during transformation
           must be calculated.
        """
        self.origin = self.idorigin.copy()

        if self.pre_align["roi"] is not None:
            self.origin[1] += self.pre_align["roi"][0][0]
            self.origin[0] += self.pre_align["roi"][1][0]

        if self.pre_transform["pad"] or self.post_transform["crop"]:
            self.origin[1] += self.extend[0][0]
            self.origin[0] += self.extend[1][0]

    def get_transformation(self):
        raise NotImplementedError()

    def gettransformation(self,i,pairwise):
        """Get transformation parameters
        """
        transfo = self.get_transformation()

        # cofs is relative to i==iref or to i==0 when pairwise
        bchange = pairwise and self.alignonraw and i!=0

        if bchange:
            self.cofs[i,...] = np.dot(self.cofs[i-1,...],transfo)
        else:
            self.cofs[i,...] = transfo

        self.update_transformation(i,bchange)

    def set_transformation(self,cof,changed):
        raise NotImplementedError()

    def update_transformation(self,i,bchange):
        """Update transformation parameters
        """
        if self.usekernel:
            self.set_transformation(self.cofs[i,...],bchange)

    def settransformidentity(self,i):
        """Make this transformation the identity
        """
        self.cofs[i,...] = self.idcof

    def minimaltransformation(self):
        """If all transformations are known, they can be reduced to minimize the difference with the original images
        """
        # TODO: do something better
        self.cofs[:,0,2] -= np.median(self.cofs[:,0,2])
        self.cofs[:,1,2] -= np.median(self.cofs[:,1,2])

    def extendfromtransformation(self):
        """If all transformations are known, padding/cropping can be calculated
        """

        self.extend = ((0,0),(0,0))
        if not self.pre_transform["pad"] and not self.post_transform["crop"]:
            return

        # Transform corners
        xy = np.empty((3,4))
        xy[0,:] = [0,self.source.imgsize[1]-1,0,self.source.imgsize[1]-1]
        xy[1,:] = [0,0,self.source.imgsize[0]-1,self.source.imgsize[0]-1]
        xy[2,:] = [1,1,1,1]
        myminmax = np.append(np.min(xy,axis=1),np.max(xy,axis=1))
        maskmin = np.empty((3,4),dtype=bool)
        maskmin[0,:] = [True,False,True,False]
        maskmin[1,:] = [True,True,False,False]
        maskmin[2,:] = [False,False,False,False]
        maskmax = np.logical_not(maskmin)
        maskmax[2,:] = [False,False,False,False]

        for i in range(self.source.nimages):
            xy2 = np.dot(self.cofwithorigin(self.cofs[i,...]),xy)
            xy2[0,:] /= xy2[2,:]
            xy2[1,:] /= xy2[2,:]
            myminmax[0:3] = np.minimum(myminmax[0:3],np.min(xy2,axis=1))
            myminmax[3:] = np.maximum(myminmax[3:],np.max(xy2,axis=1))
        xmin = int(np.floor(myminmax[0]))
        ymin = int(np.floor(myminmax[1]))
        xmax = int(np.ceil(myminmax[3]))
        ymax = int(np.ceil(myminmax[4]))

        # Padding/cropping
        o1min = ymin
        o1max = ymax-self.source.imgsize[0]+1
        o2min = xmin
        o2max = xmax-self.source.imgsize[1]+1

        # Translation only:
        #o1min = np.floor(np.min(self.cofs[:,1,2])).astype(np.int)
        #o2min = np.floor(np.min(self.cofs[:,0,2])).astype(np.int)
        #o1max = np.ceil(np.max(self.cofs[:,1,2])).astype(np.int)
        #o2max = np.ceil(np.max(self.cofs[:,0,2])).astype(np.int)

        if self.pre_transform["pad"]:
            self.extend = ((o1max,-o1min),(o2max,-o2min))
        else:
            self.extend = ((o1min,-o1max),(o2min,-o2max))

        self.pre_transform["pad"] = np.any(np.asarray(self.extend)>0)
        self.post_transform["crop"] = np.any(np.asarray(self.extend)<0)

    def parsetransformation_beforeapplication(self):
        """Adapt transformations before applying them
        """
        self.minimaltransformation()
        self.extendfromtransformation()

    def getaxesaftertransformation(self,axes):
        """Image X and Y axes after transformation
        """
        if not self.pre_transform["pad"] and not self.post_transform["crop"]:
            return

        if self.source.stackdim==2:
            ind = [0,1]
        elif self.source.stackdim==1:
            ind = [0,2]
        else:
            ind = [1,2]

        for i in range(len(ind)):
            j = ind[i]

            nleft = self.extend[i][0]
            nright = self.extend[i][1]
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
        
        if self.pre_transform["pad"] or self.post_transform["crop"]:
            imgsize = (imgsize[0] + self.extend[0][0] + self.extend[0][1],
                       imgsize[1] + self.extend[1][0] + self.extend[1][1])
            
        self.dest.prepare(nimages,imgsize,self.dtype)

    def set_reference(self,img,previous=False):
        raise NotImplementedError()

    def doalign(self,refdatasetindex,refimageindex=None,aligntype=alignType.full):
        """Align datasets and save the result.
           Calculate transformation (alignment):
            fixed image: raw -> align prep
            moving image: raw -> align prep -> transformed
           Apply transformation:
            aligned image: raw -> transform prep -> transformed -> transform post
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
                    usealignedimage = j==refdatasetindex and \
                                      aligntype!=alignType.usetransfo and \
                                      self.nopre_align() and \
                                      self.nopre_transform() and \
                                      self.nopost_transform()
                    if usealignedimage:
                        img = imgaligned
                    else:
                        img = self.readimgraw(j,i)
                        img = self.transform(img,i)

                    self.writeimg(img,j,i)

    def align(self,refdatasetindex,refimageindex = None,onraw = False,pad = True,crop = False,redo = False,roi = None):
        """Alignment function that needs to be called 
        
        Args:
            refdatasetindex(int): stack to be used for alignment
            refimageindex(Optional(int)): fixed reference to align on (pairwise alignment when None)
            onraw(Optional(bool)): when doing pairwise alignment, use the previous raw or aligned images to align the next one on
            pad(Optional(bool)): make sure nothing is transformed outside the field of view (has priority over crop)
            crop(Optional(bool)): make sure no missing data is added
            redo(Optional(bool)): apply transformations without recalculating them
            roi(Optional(array-like)): use ony part of the image to align
        """

        pairwise = refimageindex is None
        if pairwise:
            self.alignonraw = onraw
        else:
            self.alignonraw = True

        self.pre_align["roi"] = roi
        self.pre_transform["pad"] = pad
        self.post_transform["crop"] = crop
        self.calculateorigin()

        if redo:
            self.doalign(refdatasetindex,aligntype=alignType.usetransfo)
        else:
            if pad or crop:
                self.doalign(refdatasetindex,refimageindex=refimageindex,aligntype=alignType.calctransfo)
                self.parsetransformation_beforeapplication()
                self.doalign(refdatasetindex,refimageindex=refimageindex,aligntype=alignType.usetransfo)
            else:
                self.doalign(refdatasetindex,refimageindex=refimageindex)

