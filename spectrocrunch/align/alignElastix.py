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
import SimpleITK as sitk
import numpy as np
from spectrocrunch.common.stdout import stdout_redirect

class alignElastix(align):

    def __init__(self,*args,**kwargs):
        super(alignElastix,self).__init__(*args,**kwargs)

        # Prepare alignment kernel
        self.elastix = sitk.SimpleElastix()
        #self.elastix.LogToFolder("")
        #self.elastix.LogToFolderOff()
        self.elastix.LogToConsoleOff()

        parameterMap = sitk.GetDefaultParameterMap('translation')
        self.defaultvalue = -999 # Elastix cannot handle NaN's!!!
        parameterMap["DefaultPixelValue"] = (str(self.defaultvalue),)
        self.elastix.SetParameterMap(parameterMap)

        self.fixed = None
        self.moving = None

        # Prepare transformation kernel
        self.transformix = sitk.SimpleTransformix()

    def replacecval(self,img):
        if self.cval is not np.nan or self.defaultvalue is not np.nan:
            if self.cval != self.defaultvalue:
                if self.defaultvalue is np.nan:
                    img[np.isnan(img)] = self.cval
                else:
                    img[img==self.defaultvalue] = self.cval

    def execute_transformkernel(self,img):
        """Transform image according with the transformation kernel
        """
        self.changefortransform(img.shape)
        img[np.isnan(img)] = 0
        self.transformix.SetInputImage(sitk.GetImageFromArray(img))
        with stdout_redirect():
            self.transformix.Execute()
        aligned = sitk.GetArrayFromImage(self.transformix.GetResultImage())
        self.replacecval(aligned)
        return aligned
 
    def execute_alignkernel(self,img):
        """Align image on reference
        """
        img[np.isnan(img)] = 0
        self.moving = sitk.GetImageFromArray(img)
        self.elastix.SetMovingImage(self.moving)

        try:
            with stdout_redirect():
                self.elastix.Execute()
            aligned = sitk.GetArrayFromImage(self.elastix.GetResultImage())
            self.replacecval(aligned)
        except:
            aligned = img.copy()

        return aligned

    def set_reference(self,img,previous=False):
        """Reference for alignment
        """
        if previous:
            self.fixed = self.moving    
        else:
            img[np.isnan(img)] = 0
            self.fixed = sitk.GetImageFromArray(img)
        self.elastix.SetFixedImage(self.fixed)

    def get_transformation(self):
        """Get transformation from alignment kernel.
        """
        transformParameterMap = self.elastix.GetTransformParameterMap()
        # Elastix takes x as the first dimension while Python needs y as the first dimension.
        return np.array(transformParameterMap[0]["TransformParameters"][::-1], self.dtype)

    def set_transformation(self,offset,changed):
        """Set the transformation kernel according to the alignment kernel and adapted transformation
        """
        transformParameterMap = self.elastix.GetTransformParameterMap()
        if changed:
            # Elastix takes x as the first dimension while Python needs y as the first dimension.
            transformParameterMap[0]["TransformParameters"] = (str(offset[1]),str(offset[0]))
        self.transformix.SetTransformParameterMap(transformParameterMap)

    def changefortransform(self,shape):
        transformParameterMap = self.elastix.GetTransformParameterMap()

        oldshape = (int(transformParameterMap[0]["Size"][1]),int(transformParameterMap[0]["Size"][0]))
        if shape==oldshape:
            return

        # Not sure about this
        transformParameterMap[0]["Size"] = (str(shape[1]),str(shape[0]))
        transformParameterMap[0]["Origin"] = (str(self.origin[1]),str(self.origin[0]))

        self.transformix.SetTransformParameterMap(transformParameterMap)


