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


from .align import align
from .types import transformationType
import numpy as np
import scipy.ndimage

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
            xy = scipy.ndimage.measurements.center_of_mass(self.handle_missing(img,0))
        elif self.xytype=="min":
            xy = np.unravel_index(np.nanargmin(self.handle_missing(img,np.nan)),img.shape)
        else: # self.xytype=="max"
            xy = np.unravel_index(np.nanargmax(self.handle_missing(img,np.nan)),img.shape)
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
            self._transform.set(transform)

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

