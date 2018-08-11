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

import numpy as np
from matplotlib import pyplot as plt
from shapely.geometry import Polygon
from shapely.geometry import LineString

from .alignSource import alignSource
from .alignDest import alignDest
from .types import alignType
from .types import transformationType
from .transform import transform
from ..common.roi import cliproi

import logging

class align(object):
    """Allows for the alignment of several stacks based on one stack.
       "Alignment" is the process of determining a transformation between
       two images that should represent the same thing but transformed/deformed.

        Transformation:
            L.Xold = Xnew (coordinate transformation)
            Xold = C.Xnew (change-of-frame)
        Combine transformations:
            L2.L1.Xold = Xnew (coordinate transformation, multiply from the left)
            Xold = C1.C2.Xnew (change-of-frame, multiply from the right)
        Effect of change-of-frame:
            L2 = C^-1.L1.C (coordinate transformation in new frame)
            C2 = C^-1.C1.C (change-of-frame in new frame)

        Calculate transformation (image alignment):
            fixed image:   raw -> dopre_align (cof=C1) -> img1 
            moving image:  raw -> dopre_align (cof=C1) -> img2 -> execute_alignkernel(img1,img2) -> img3, C2

            C1^-1: pre-align to raw (roi)
            C2: alignment in pre-align frame
            C1.C2.C1^-1: cof in raw frame
                           
        Apply transformation (image transformation):
            aligned image: raw -> dopre_transform (cof=C3) -> img1 -> execute_transform() -> img4 -> dopost_transform (cof=C4) -> img5

            C3: raw to pre-transform (pad)
            C1^-1.C3: pre-align to pre-transform
            C3^-1.C1.C2.C1^-1.C3: alignment in pre-transform frame
            C4: post transformation (crop)
    """

    def __init__(self,source,sourcelist,dest,destlist,extension,\
                stackdim=None,overwrite=False,cval=np.nan,plot=False,\
                transfotype=transformationType.translation):
        """
        """

        # Data IO
        self.source = alignSource(source,sourcelist,stackdim=stackdim)
        self.dest = alignDest(dest,destlist,extension,stackdim=self.source.stackdim,overwrite=overwrite)
    
        # Missing data
        self.cval = cval

        # Transformation settings (set before actual transformation)
        self.dtype = (np.float32(1)*self.source.dtype.type(1)).dtype.type
        self.alignonraw = True
        self.usekernel = False # Doesn't work well for Elastix!
        self.pre_align = {"roi":None}
        self.pre_transform = {"pad":False}
        self.post_transform = {"crop":False}
        self.pre_transform_requested = {"pad":False}
        self.post_transform_requested = {"crop":False}
        self.extend = ((0,0),(0,0)) # negative: crop, positive: pad

        # Transformation (change of frame matrices, not change of coordinates!)
        self.transfotype = transfotype
        self.transfos = [self.defaulttransform() for i in range(self.source.nimages)]
        self.C1 = self.defaulttransform(ttype=transformationType.translation)
        self.C1inv = self.defaulttransform(ttype=transformationType.translation)
        self.C3invC1 = self.defaulttransform(ttype=transformationType.translation)
        self.C1invC3 = self.defaulttransform(ttype=transformationType.translation)

        # Plot
        self.plotinfo = {"ON":plot,"fig":None,"axes":None}

    def defaulttransform(self,ttype=None,dtype=None):
        if ttype is None:
            ttype = self.transfotype
        if dtype is None:
            dtype = self.dtype
        return transform(ttype,dtype=dtype,cval=self.cval)

    def enableplot(self):
        self.plotinfo["ON"] = True

    def disableplot(self):
        self.plotinfo["ON"] = False

    def plot(self,img,index,title):
        """Visualize alignment in progress
        """
        if not self.plotinfo["ON"]:
            return

        if self.plotinfo["fig"] is None:
            self.plotinfo["fig"],self.plotinfo["axes"] = plt.subplots(1,3)
        ax = self.plotinfo["axes"][index]
        ax.cla()
        
        if img.size in img.shape:
            ax.plot(img.flatten())
        else:
            #img2 = img.copy()
            #img2[np.isnan(img2)] = 0
            ax.imshow(img,origin='lower',interpolation='nearest',cmap='jet')
        ax.set_title(title)

        plt.pause(0.01)
        
    def padfromextend(self):
        return ((max(self.extend[0][0],0),max(self.extend[0][1],0)),\
               (max(self.extend[1][0],0),max(self.extend[1][1],0)))
               
    def cropfromextend(self,dim1,dim2):
        return  ((max(-self.extend[0][0],0),dim1-max(-self.extend[0][1],0)),\
                (max(-self.extend[1][0],0),dim2-max(-self.extend[1][1],0)))
                
    def pad(self,img):
        """Apply padding
        """
        pad = self.padfromextend()
        if np.count_nonzero(pad)!=0:
            return np.pad(img,pad,'constant',constant_values=(self.cval,self.cval))
        else:
            return img

    def crop(self,img):
        """Apply cropping
        """
        dim1,dim2 = img.shape
        crop = self.cropfromextend(dim1,dim2)
        if crop[0][0]!=0 or crop[1][0]!=0 or crop[0][1]!=dim1 or crop[1][1]!=dim2:
            return img[crop[0][0]:crop[0][1],crop[1][0]:crop[1][1]]
        else:
            return img

    def roi(self,img,roi):
        """Extract ROI
        """
        [[ya,yb],[xa,xb]] = cliproi(img.shape,roi)
        if xb<=xa or yb<=ya:
            raise ValueError("ROI reduces image size to zero: [{}:{},{}:{}]".format(ya,yb,xa,xb))
        return img[ya:yb,xa:xb]

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
        if 0 in img.shape or len(img.shape)!=2:
            raise ValueError("Image preprocessed for alignment has shape {}".format(img.shape))
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
        if self.pre_transform["pad"]:
            img = self.pad(img)
        return img

    def nopost_transform(self):
        """
        Returns:
            bool: True when transformation doesn't have any post processing
        """
        return not self.post_transform["crop"]

    def dopost_transform(self,img):
        """Manual transformation after the real transformation (not used in alignment)
        """
        if self.post_transform["crop"]:
            img = self.crop(img)
        return img

    def execute_alignkernel(self,img):
        raise NotImplementedError()

    def execute_transformkernel(self,img):
        raise NotImplementedError()

    def execute_transform_nokernel(self,img,C2):
        """Apply a transformation to an image
        """
        return self.cof_in_pretransform_frame(C2).transformimage(img)

    def absolute_cofs(self,homography=False):
        """Change-of-frame for each raw image (i.e. maps destination keypoints to reference keypoints)
        """
        if self.source.nimages==0:
            return None
        if homography:
            return np.array([t.getnumpyhomography() for t in self.transfos])
        else:
            return np.array([t.getnumpy() for t in self.transfos])

    def cof_in_raw_frame(self,C2):
        """ C2' = Ch^-1.C2.Ch
            C2: in pre-align frame
            C2': in raw frame
            Ch: pre-align frame to raw frame
        """
        if self.C1inv.isidentity():
            return C2
        return self.C1.dot(C2).dot(self.C1inv)

    def cof_in_pretransform_frame(self,C2):
        """ C2' = Ch^-1.C2.Ch
            C2: in pre-align frame
            C2': in pre-transform frame
            Ch: pre-align frame to pre-transform frame
        """
        if self.C1invC3.isidentity():
            return C2
        return self.C3invC1.dot(C2).dot(self.C1invC3)

    def calccof_raw_to_prealign(self):
        """ raw -> pre-align: C1 (roi)
            pre-align -> raw: C1^-1
        """
        self.C1.setidentity()
        self.C1inv.setidentity()
        if self.pre_align["roi"] is not None:
            self.C1.settranslation([self.pre_align["roi"][1][0],self.pre_align["roi"][0][0]])
            self.C1inv.settranslation([-self.pre_align["roi"][1][0],-self.pre_align["roi"][0][0]])
        self.C3invC1.fromtransform(self.C1)
        self.C1invC3.fromtransform(self.C1inv)

    def calccof_pretransform_to_prealign(self):
        """ raw -> pre-align: C1 (roi)
            raw -> pre-transform: C3 (pad)
            pre-transform -> pre-align: C3^-1.C1
            pre-align -> pre-transform: C1^-1.C3
        """
        self.C3invC1.fromtransform(self.C1)
        self.C1invC3.fromtransform(self.C1inv)

        if self.pre_transform["pad"]:
            C3inv = self.defaulttransform(ttype=transformationType.translation)
            C3inv.settranslation([max(self.extend[1][0],0),max(self.extend[0][0],0)]) # o2min, o1min
            C3 = C3inv.inverse()
            self.C3invC1.dotleftinplace(C3inv)
            self.C1invC3.dotinplace(C3)

            #TODO: not this right?
            #self.C3invC1.settranslation([self.extend[1][0],self.extend[0][0]])
            #self.C1invC3.settranslation([-self.extend[1][0],-self.extend[0][0]])

    def execute_transform(self,img,i):
        """Transform according to the transformation extracted from the transformation kernel (see gettransformation).
        """
        if not self.transformidentity(i):
            if self.usekernel:
                return self.execute_transformkernel(img)
            else:
                return self.execute_transform_nokernel(img,self.transfos[i])
        return img

    def transformidentity(self,i):
        """Is the transformation the identity
        """
        return self.transfos[i].isidentity()

    def pureidentity(self,i):
        """Is the transformation the identity, including the changes applied before (padding) and after (cropping)
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
        imgtransformed = self.execute_transform(imgtransformed,i)

        #plt.figure(5)
        #plt.imshow(imgtransformed,origin='lower',interpolation='nearest',cmap='jet')
        #plt.pause(1)

        # Apply final transformation (not used in alignment)
        imgtransformed = self.dopost_transform(imgtransformed)

        return imgtransformed

    def get_transformation(self):
        """To be implemented by a derived class
        """
        raise NotImplementedError()

    def gettransformation(self,i,pairwise):
        """Get transformation parameters
        """
        transfo = self.get_transformation()

        # cofs are relative to i==iref or to i==0 when pairwise
        bchange = pairwise and self.alignonraw and i!=0

        if bchange:
            self.transfos[i].fromtransform(self.transfos[i-1].dot(transfo)) # new cof multiplied at the right
        else:
            self.transfos[i].fromtransform(transfo)

        self.update_transformation(i,bchange)

    def set_transformation(self,cof,changed):
        """To be implemented by a derived class
        """
        raise NotImplementedError()

    def update_transformation(self,i,bchange):
        """Set the active transformation
        """
        if self.usekernel:
            if self.C1invC3.isidentity():
                self.set_transformation(self.transfos[i],bchange)
            else:
                self.set_transformation(self.cof_in_pretransform_frame(self.transfos[i]),True)

    def settransformidentity(self,i):
        """Make this transformation the identity
        """
        self.transfos[i].setidentity()

    def genpolygon(self,lst):
        p = Polygon(lst)
        if p.area==0:
            p = LineString(lst)
        return p

    def polygonempty(self,p):
        if isinstance(p,Polygon):
            return p.area==0
        else:
            return p.length==0

    def untransformedimagepolygon(self):
        add0 = 0.
        add1 = 0.
        return self.genpolygon([(add0,add0),\
                        (self.source.imgsize[1]-1+add1,add0),\
                        (self.source.imgsize[1]-1+add1,self.source.imgsize[0]-1+add1),\
                        (add0,self.source.imgsize[0]-1+add1)])

    def transformedimagepolygons(self):
        add0 = 0.
        add1 = 0.

        # Corners of the image in the transformed frame: A'
        xy = np.empty((3,4))
        xy[0,:] = [add0,self.source.imgsize[1]-1+add1,self.source.imgsize[1]-1+add1,add0] #x
        xy[1,:] = [add0,add0,self.source.imgsize[0]-1+add1,self.source.imgsize[0]-1+add1] #y
        xy[2,:] = [1,1,1,1]

        # Corners of the image in the raw frame: A = C.A'
        ret = [None]*self.source.nimages
        for i in range(self.source.nimages):
            xy2 = self.cof_in_raw_frame(self.transfos[i]).transformcoordinates(xy) # C1^(-1).C2^(-1).C1.XY
            xy2[0,:] /= xy2[2,:]
            xy2[1,:] /= xy2[2,:]
            ret[i] = self.genpolygon(xy2[0:2,:].T)
    
        return ret

    def polygoncollectionbounds(self,ps,pad=True):
        p = ps[0]
        if pad:
            for i in range(1,len(ps)):
                p = p.union(ps[i])
        else:
            for i in range(1,len(ps)):
                p = p.intersection(ps[i])

            if self.polygonempty(p):
                logger = logging.getLogger(__name__)
                logger.warning("Cropping skipped because there would be nothing left.")
                return ()

        xmin, ymin, xmax, ymax = p.bounds
        xmin = int(np.floor(xmin))
        ymin = int(np.floor(ymin))
        xmax = int(np.ceil(xmax))
        ymax = int(np.ceil(ymax))

        return xmin, ymin, xmax, ymax

    def minimaltransformation(self,p0,ps,centroids=False):
        """If all transformations are known, they can be reduced to minimize the difference with the original image
        """
        #return

        if centroids:
            # Put average centroids in the middle
            x0,y0 = p0.centroid.coords.xy
            xy = np.array([p.centroid.coords.xy for p in ps])
            x0 = x0[0]
            y0 = y0[0]
            x = np.mean(xy[:,0])
            y = np.mean(xy[:,1])
        else:
            # Center total boundary
            xmin0, ymin0, xmax0, ymax0 = p0.bounds
            xmin, ymin, xmax, ymax = self.polygoncollectionbounds(ps)
            x0 = np.mean([xmin0,xmax0])
            y0 = np.mean([ymin0,ymax0])
            x = np.mean([xmin,xmax])
            y = np.mean([ymin,ymax])
        
        # Center
        trn = self.defaulttransform(ttype=transformationType.translation)
        trn.settranslation(x-x0,y-y0)
        for t in self.transfos:
            t.dotinplace(trn)
            
    def setextend(self,xmin, ymin, xmax, ymax):
        # Padding/cropping <> positive/negative
        o1min = -ymin
        o1max = ymax-self.source.imgsize[0]+1
        o2min = -xmin
        o2max = xmax-self.source.imgsize[1]+1
        
        self.extend = ((o1min,o1max),(o2min,o2max))
        
        self.pre_transform["pad"] = np.any(np.asarray(self.extend)>0)
        self.post_transform["crop"] = np.any(np.asarray(self.extend)<0)
        
    def extendfromtransformation(self,p0,ps):
        """If all transformations are known, padding/cropping can be calculated
        """
        self.extend = ((0,0),(0,0))

        # Smallest rectangle that contains the union (pad) or intersection (crop) of all polygons
        tmp = self.polygoncollectionbounds(ps,pad=self.pre_transform_requested["pad"])
        if len(tmp)!=4:
            self.pre_transform["pad"] = False
            self.post_transform["crop"] = False
            return
            
        self.setextend(*tmp)
        
    def bextendfrommask(self,pairwise):
        return self.post_transform_requested["crop"] and\
               self.cval is np.nan and\
               (self.pre_align["roi"] is None or self.transfotype==transformationType.translation) and\
               not pairwise
        
    def setextendmask(self,img,reset=False):
        mask = np.logical_not(np.isnan(img))   
        #if self.cval is np.nan:
        #    mask = np.logical_not(np.isnan(img))
        #else:
        #    mask = img != self.cval
        
        if reset:
            self._extendmask = mask
        else:
            self._extendmask &= mask

    def extendfrommask(self):
        """If all transformations are applied, padding/cropping can be calculated
        """
        indvalidrow = np.argwhere(self._extendmask.sum(axis=1))
        indvalidcol = np.argwhere(self._extendmask.sum(axis=0))
        
        # When pre_align["roi"]: only valid for translations
        ymin = indvalidrow[0][0]
        ymax = indvalidrow[-1][0]-self._extendmask.shape[0]+self.source.imgsize[0]
        xmin = indvalidcol[0][0]
        xmax = indvalidcol[-1][0]-self._extendmask.shape[1]+self.source.imgsize[1]
        
        self.setextend(xmin,ymin,xmax,ymax)
        
    def parsetransformation_beforeapplication(self,pairwise):
        """Adapt transformations before applying them
        """
        if self.bextendfrommask(pairwise):
            self.extendfrommask()
        else:
            # Corners of the image in the transformed frame
            p0 = self.untransformedimagepolygon()

            # Corners of the transformed image in the raw frame
            ps = self.transformedimagepolygons()

            # Adapt transformation
            if self.pre_transform_requested["pad"] or self.post_transform_requested["crop"]:
                # adapt self.extend to either crop or pad
                self.extendfromtransformation(p0,ps)
            else:
                # try to fit as much data in the original image size as possible
                self.minimaltransformation(p0,ps)
                
        self.calccof_pretransform_to_prealign()

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

    def dest_imgsize(self,nopost=False):
        imgsize = self.source.imgsize
        if self.pre_transform["pad"] or (self.post_transform["crop"] and not nopost):
            imgsize = (imgsize[0] + self.extend[0][0] + self.extend[0][1],
                       imgsize[1] + self.extend[1][0] + self.extend[1][1])
        return imgsize
        
    def preparedestination(self,img=None):
        """Allocate space for saving results
        """
        if img is not None:
            self.setup_post_transform(img)
        
        nimages = self.source.nimages
        imgsize = self.dest_imgsize()
        self.dest.prepare(nimages,imgsize,self.dtype)

    def set_reference(self,img,previous=False):
        raise NotImplementedError()

    def doalign(self,refdatasetindex,refimageindex=None,aligntype=alignType.full):
        """Align datasets and save the result.
        """
        logger = logging.getLogger(__name__)
 
        pairwise = refimageindex is None

        # Prepare destination (will be done after alignment)
        if aligntype!=alignType.calctransfo:
            self.preparedestination()

        # First reference image
        if aligntype!=alignType.usetransfo:
            if pairwise:
                # Pair-wise alignment: first image is the first reference
                imgref = self.readimgrawprep(refdatasetindex,0)
                iref = 0
            else:
                # Fixed-reference alignment
                rawprep = self.readimgrawprep(refdatasetindex,refimageindex)
                iref = refimageindex
                self.plot(rawprep,0,"Reference %d (fixed)"%iref)
                self.set_reference(rawprep)

        #from pympler import tracker
        #tr = tracker.SummaryTracker()
        #s1 = None

        # Loop over the images
        bfirst = True
        for i in range(self.source.nimages):
            if aligntype!=alignType.usetransfo:
                # Image i
                rawprep = self.readimgrawprep(refdatasetindex,i)
                #np.save("img{}.npy".format(i),rawprep)
                
                # Update fixed image
                if pairwise:
                    self.set_reference(imgref)
                
                # Get align transformation
                logger.debug("Align index %d on %d"%(i,iref))
                
                if i == iref:
                    self.settransformidentity(i)
                    imgaligned = rawprep
                else:
                    # Align image i to reference

                    #if s1 is None:
                    #    s1 = tr.create_summary()

                    imgaligned = self.execute_alignkernel(rawprep)
                    #s2 = tr.create_summary()

                    #tr.print_diff(summary1=s1,summary2=s2)
                    if pairwise:
                        self.plot(imgref,0,"Reference %d (pair-wise)"%(iref))
                    self.plot(rawprep,2,"To align %d"%i)
                    self.plot(imgaligned,1,"Aligned %d"%i)
                    self.gettransformation(i,pairwise)
                
                logger.debug("Aligned index %d on %d"%(i,iref))

                # Reference for the next image
                if pairwise:
                    if self.alignonraw:
                        imgref = rawprep
                    else:
                        imgref = imgaligned
                    iref = i

                # Only calculation
                if aligntype==alignType.calctransfo:
                    # This is still not good enough because different kernels are used
                    # when calculating anf applying alignment
                    if self.bextendfrommask(pairwise):
                        self.setextendmask(imgaligned,reset=bfirst)
                    bfirst = False
                    continue # no results needed
                        
            # Save the transformed image i of all datasets
            if self.pureidentity(i):
                for j in range(self.source.nsets):
                    self.copyimg(j,i)
            else:
                if aligntype==alignType.usetransfo:
                    self.update_transformation(i,True)

                for j in range(self.source.nsets):
                    usealignedimage = aligntype!=alignType.usetransfo and \
                                      j==refdatasetindex and \
                                      self.nopre_align() and \
                                      self.nopre_transform() and \
                                      self.nopost_transform()
                    if usealignedimage:
                        img = imgaligned
                    else:
                        img = self.readimgraw(j,i)
                        img = self.transform(img,i)

                    self.writeimg(img,j,i)
            
    def setroi(self,roi):
        if roi is None:
            self.pre_align["roi"] = None
        else:
            self.pre_align["roi"] = ((0 if roi[0][0] is None else roi[0][0],roi[0][1]),\
                                     (0 if roi[1][0] is None else roi[1][0],roi[1][1]))
        self.calccof_raw_to_prealign()

    def align(self,refdatasetindex,refimageindex = None,onraw = False,pad = True,crop = False,redo = False,roi = None):
        """Alignment function that needs to be called 
        
        Args:
            refdatasetindex(int): stack to be used for alignment
            refimageindex(Optional(int)): fixed reference to align on (pairwise alignment when None)
            onraw(Optional(bool)): when doing pairwise alignment, use the previous raw or aligned images to align the next one on
            pad(Optional(bool)): make sure nothing is transformed outside the field of view (has priority over crop)
            crop(Optional(bool)): make sure no missing data is added
            redo(Optional(bool)): apply transformations without recalculating them (all other keywords are ignored)
            roi(Optional(array-like)): use ony part of the image to align
        """
        if redo:
            self.doalign(refdatasetindex,aligntype=alignType.usetransfo)
        else:
            pairwise = refimageindex is None
            if pairwise:
                self.alignonraw = onraw
            else:
                self.alignonraw = True

            self.setroi(roi)
            self.pre_transform_requested["pad"] = pad
            self.post_transform_requested["crop"] = crop
            self.pre_transform["pad"] = pad
            self.post_transform["crop"] = crop
            center = False
        
            if roi or pad or crop or center:
                self.doalign(refdatasetindex,refimageindex=refimageindex,aligntype=alignType.calctransfo)
                self.parsetransformation_beforeapplication(pairwise)
                self.doalign(refdatasetindex,refimageindex=refimageindex,aligntype=alignType.usetransfo)
            else:
                self.doalign(refdatasetindex,refimageindex=refimageindex)

