# -*- coding: utf-8 -*-

from .align import align
from .types import transformationType
import numpy as np
from ..math import center

class alignSimple(align):

    def __init__(self,*args,**kwargs):
        super(alignSimple,self).__init__(*args,**kwargs)

        # Images
        self.fixedxy = None
        self.movingxy = None

        # change of reference frame
        self._transform = self.defaulttransform()

    def execute_transformkernel(self,img):
        """Transform image according with the transformation kernel
        """
        return self.execute_transform_nokernel(img,self._transform)

    def execute_alignkernel(self,img):
        """Align image on reference
        """
        if self.transfotype!=transformationType.translation:
            raise NotImplementedError("Sift doesn't support this type of transformation.")

        self.movingxy = self.getxy(img)
        self._transform.settranslation(self.movingxy-self.fixedxy)
        return self.execute_transformkernel(img)

    def handle_missing(self,img,newval):
        """Handling of missing data
        """
        if self.cval is np.nan and newval is np.nan:
            return img
        if self.cval == newval:
            return img

        if self.cval != 0:
            if self.cval is np.nan:
                missing = np.isnan(img)
            else:
                missing = img==self.cval
            bmissing = np.any(missing)
        else:
            bmissing = False

        if bmissing:
            img2 = img.copy()
            img2[missing] = newval
        else:
            img2 = img

        return img2

    def getxy(self,img):
        """Get marker (min, max, centroid)
        """
        xy = None
        if self.xytype=="centroid":
            xy = center.fcentroid(self.handle_missing(img,0))
        elif self.xytype=="min":
            xy = center.fmin(self.handle_missing(img,np.nan))
        elif self.xytype=="gaussmax":
            xy = center.fgaussmax(self.handle_missing(img,np.nan))
        else: # self.xytype=="max"
            xy = center.fmax(self.handle_missing(img,np.nan))

        if img.size in img.shape:
            if img.shape[0]==1:
                xy = (0,xy)
            else:
                xy = (xy,0)

        return np.array(xy)[::-1]

    def set_reference(self,img,previous=False):
        """Reference for alignment
        """
        if previous:
            self.fixedxy = self.movingxy
        else:
            self.fixedxy = self.getxy(img)

    def get_transformation(self):
        """Get transformation
        """
        return self._transform

    def set_transformation(self,cof,changed):
        """Set transformation
        """
        if changed:
            self._transform.fromtransform(transform)

class alignMin(alignSimple):
    def __init__(self,*args,**kwargs):
        super(alignMin,self).__init__(*args,**kwargs)
        self.xytype = "min"

class alignMax(alignSimple):
    def __init__(self,*args,**kwargs):
        super(alignMax,self).__init__(*args,**kwargs)
        self.xytype = "max"

class alignCentroid(alignSimple):
    def __init__(self,*args,**kwargs):
        super(alignCentroid,self).__init__(*args,**kwargs)
        self.xytype = "centroid"

class alignGaussMax(alignSimple):
    def __init__(self,*args,**kwargs):
        super(alignGaussMax,self).__init__(*args,**kwargs)
        self.xytype = "gaussmax"

