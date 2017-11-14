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

from .types import transformationType
import numpy as np
import scipy.ndimage.interpolation

class transform(object):

    def __init__(self,transfotype,dtype=np.float,cval=np.nan,**interpolationargs):
        self.transfotype = transfotype
        self.dtype = dtype
        self.cval = cval

        if "order" not in interpolationargs:
            interpolationargs["order"] = 1
        if "mode" not in interpolationargs:
            interpolationargs["mode"] = "constant"
        self.interpolationargs = interpolationargs
        
        # (change of frame matrices, not change of coordinates!)
        self.cof = self.getidentity()

    def set(self,transformobj):
        self.transfotype = transformobj.transfotype
        self.dtype = transformobj.dtype
        self.cval = transformobj.cval
        self.cof[:] = transformobj.cof[:]

    def getnumpy(self):
        return self.cof[:]

    def getnumpyhomography(self):
        if self.transfotype==transformationType.translation:
            return np.array([[1,0,self.cof[0]],[0,1,self.cof[1]],[0,0,1]],dtype=self.dtype)
        else:
            return self.cof[:]

    def settranslation(self,trn,ty=None):
        """
        Args:
            trn(number or 2-array-like): x and optionally y translation
            ty(Optional(number)): y translation
        """
        if ty is not None:
            trn = [trn,ty]
        if self.transfotype==transformationType.translation:
            self.cof[:] = trn
        elif self.transfotype==transformationType.rigid:
            self.cof[0:2,2] = trn
        elif self.transfotype==transformationType.similarity:
            self.cof[0:2,2] = trn
        elif self.transfotype==transformationType.affine:
            self.cof[0:2,2] = trn
        elif self.transfotype==transformationType.homography:
            self.cof[0:2,2] = trn
        else:
            raise ValueError("Transformation does not have an translation part")

    def gettranslation(self):
        if self.transfotype==transformationType.translation:
            return self.cof[:]
        elif self.transfotype==transformationType.rigid:
            return self.cof[0:2,2]
        elif self.transfotype==transformationType.similarity:
            return self.cof[0:2,2]
        elif self.transfotype==transformationType.affine:
            return self.cof[0:2,2]
        elif self.transfotype==transformationType.homography:
            return self.cof[0:2,2]
        else:
            raise ValueError("Transformation does not have an translation part")

    def setrigid(self,theta,tx,ty):
        if self.transfotype==transformationType.translation:
            trn = self.gettranslation()
            self.transfotype==transformationType.affine            
            self.cof = self.getidentity()
        elif self.transfotype==transformationType.rigid:
            pass
        elif self.transfotype==transformationType.similarity:
            self.transfotype=transformationType.rigid
        elif self.transfotype==transformationType.affine:
            self.transfotype=transformationType.rigid
        elif self.transfotype==transformationType.homography:
            pass
        else:
            raise ValueError("Transformation does not have an linear part")

        costheta = np.cos(theta)
        sintheta = np.sin(theta)
        self.cof[0:2,:] = np.array([[costheta,-sintheta,tx],[sintheta,costheta,ty]], self.dtype)

    def getrigid(self):
        if self.transfotype==transformationType.translation:
            return 0,self.cof[0],self.cof[1]
        elif self.transfotype==transformationType.rigid:
            return np.arccos(self.cof[0,0]),self.cof[0,2],self.cof[1,2]
        else:
            raise ValueError("Transformation should be a translation or a rigid transformation")

    def setlinear(self,M):
        if self.transfotype==transformationType.translation:
            trn = self.gettranslation()
            self.transfotype=transformationType.affine            
            self.cof = self.getidentity()
            self.settranslation(trn)
        elif self.transfotype==transformationType.rigid:
            pass
        elif self.transfotype==transformationType.similarity:
            pass
        elif self.transfotype==transformationType.affine:
            pass
        elif self.transfotype==transformationType.homography:
            pass
        else:
            raise ValueError("Transformation does not have an linear part")

        self.cof[0:2,0:2] = M

    def getlinear(self):
        if self.transfotype==transformationType.translation:
            return np.identity(2,dtype = self.dtype)
        elif self.transfotype==transformationType.rigid:
            return self.cof[0:2,0:2]
        elif self.transfotype==transformationType.similarity:
            return self.cof[0:2,0:2]
        elif self.transfotype==transformationType.affine:
            return self.cof[0:2,0:2]
        elif self.transfotype==transformationType.homography:
            return self.cof[0:2,0:2]
        else:
            raise ValueError("Transformation does not have an linear part")

    def setaffine(self,M):
        if self.transfotype==transformationType.translation:
            trn = self.gettranslation()
            self.transfotype=transformationType.affine            
            self.cof = self.getidentity()
            self.settranslation(trn)
        elif self.transfotype==transformationType.rigid:
            pass
        elif self.transfotype==transformationType.similarity:
            pass
        elif self.transfotype==transformationType.affine:
            pass
        elif self.transfotype==transformationType.homography:
            pass
        else:
            raise ValueError("Transformation does not have an affine part")

        self.cof[0:2,:] = M

    def getaffine(self):
        if self.transfotype==transformationType.translation:
            return np.array([[1,0,self.cof[0]],[0,1,self.cof[1]]],dtype=self.dtype)
        elif self.transfotype==transformationType.rigid:
            return self.cof[:,0:2]
        elif self.transfotype==transformationType.similarity:
            return self.cof[:,0:2]
        elif self.transfotype==transformationType.affine:
            return self.cof[:,0:2]
        elif self.transfotype==transformationType.homography:
            return self.cof[:,0:2]
        else:
            raise ValueError("Transformation does not have an affine part")

    def sethomography(self,M):
        if self.transfotype==transformationType.translation:
            trn = self.gettranslation()
            self.transfotype=transformationType.homography            
            self.cof = self.getidentity()
            self.settranslation(trn)
        elif self.transfotype==transformationType.rigid:
            self.transfotype==transformationType.homography
        elif self.transfotype==transformationType.similarity:
            self.transfotype==transformationType.homography
        elif self.transfotype==transformationType.affine:
            self.transfotype==transformationType.homography
        elif self.transfotype==transformationType.homography:
            pass
        else:
            raise ValueError("Transformation does not have an affine part")

        self.cof[:] = M

    def gethomography(self):
        if self.transfotype==transformationType.translation:
            return np.array([[1,0,self.cof[0]],[0,1,self.cof[1]],[0,0,1]],dtype=self.dtype)
        elif self.transfotype==transformationType.rigid:
            return self.cof[:]
        elif self.transfotype==transformationType.similarity:
            return self.cof[:]
        elif self.transfotype==transformationType.affine:
            return self.cof[:]
        elif self.transfotype==transformationType.homography:
            return self.cof[:]
        else:
            raise ValueError("Transformation does not have an affine part")

    def setidentity(self):
        if self.transfotype==transformationType.translation:
            self.cof[:] = np.zeros(2,dtype = self.dtype)
        else:
            self.cof[:] = np.identity(3,dtype = self.dtype)

    def getidentity(self):
        if self.transfotype==transformationType.translation:
            return np.zeros(2,dtype = self.dtype)
        else:
            return np.identity(3,dtype = self.dtype)

    def isidentity(self):
        if self.transfotype==transformationType.translation:
            return np.array_equal(self.cof,np.zeros(2,dtype = self.dtype))
        else:
            return np.array_equal(self.cof,np.identity(3,dtype = self.dtype))  

    def islinearidentity(self):
        if self.transfotype==transformationType.translation:
            return True
        else:
            return np.array_equal(self.cof[0:2,0:2],np.identity(2,dtype = self.dtype))  

    def isprojidentity(self):
        if self.transfotype==transformationType.translation:
            return True
        else:
            return np.array_equal(self.cof[2,0:2],np.zeros(2,dtype = self.dtype))

    def _dot(self,C,right=True):
        if self.transfotype==transformationType.translation and C.transfotype==transformationType.translation:
            cof = self.cof + C.cof
            return cof,transformationType.translation

        if self.transfotype==transformationType.translation:
            C1 = np.identity(3,dtype = self.dtype)
            C1[0:2,2] = self.cof
        else:
            C1 = self.cof

        if C.transfotype==transformationType.translation:
            C2 = np.identity(3,dtype = self.dtype)
            C2[0:2,2] = C.cof
        else:
            C2 = C.cof

        if right:
            cof = np.dot(C1,C2)
        else:
            cof = np.dot(C2,C1)

        if self.transfotype==transformationType.homography or C.transfotype==transformationType.homography:
            transfotype=transformationType.homography
        elif self.transfotype==transformationType.affine or C.transfotype==transformationType.affine:
            transfotype=transformationType.affine
        elif self.transfotype==transformationType.similarity or C.transfotype==transformationType.similarity:
            transfotype=transformationType.similarity
        else:
            transfotype=transformationType.rigid

        return cof,transfotype

    def dot(self,C):
        cof,transfotype = self._dot(C)
        ret = transform(transfotype,dtype=cof.dtype,cval=self.cval)
        ret.cof[:] = cof
        return ret

    def dotinplace(self,C):
        # C after self
        cof,transfotype = self._dot(C)
        if len(cof)==len(self.cof):
            self.cof[:] = cof
        else:
            self.cof = cof
        self.transfotype = transfotype
        self.dtype = cof.dtype

    def dotleftinplace(self,C):
        # C before self
        cof,transfotype = self._dot(C,right=False)
        if len(cof)==len(self.cof):
            self.cof[:] = cof
        else:
            self.cof = cof
        self.transfotype = transfotype
        self.dtype = cof.dtype

    def inverse(self):
        ret = transform(self.transfotype,dtype=self.cof.dtype,cval=self.cval)
        if self.transfotype==transformationType.translation:
            ret.cof[:] = -self.cof
        else:
            ret.cof[:] = np.linalg.inv(self.cof)
        return ret

    def inverseinplace(self):
        if self.transfotype==transformationType.translation:
            self.cof[:] = -self.cof
        else:
            self.cof[:] = np.linalg.inv(self.cof)

    def transformcoordinates(self,xy):
        # Anew = C^-1.Aold
        if self.transfotype==transformationType.translation:
            Ci = -self.cof
        else:
            Ci = np.linalg.inv(self.cof)
        return self._transformcoordinates(Ci,xy)

    def transformcoordinatesinverse(self,xy):
        # Aold = C.Anew
        return self._transformcoordinates(self.cof,xy)

    def _transformcoordinates(self,C,xy):
        # xy' = C.xy
        if self.transfotype==transformationType.translation:
            if len(xy.shape)==1:
                if xy.shape[0]==2:
                    return xy + C
                elif xy.shape[0]==3:
                    return np.append(xy[0:2] + C,xy[2])
                else:
                    raise ValueError("Shape of coordinates cannot be handled")
            else:
                if xy.shape[0]==2:
                    return xy + C[:,np.newaxis]
                elif xy.shape[0]==3:
                    return np.vstack((xy[0:2,:] + C[:,np.newaxis],xy[2,:]))
                else:
                    raise ValueError("Shape of coordinates cannot be handled")
        else:
            return np.dot(C,xy)

    def transformimage(self,img):
        # self.cof: change-of-frame matrix for coordinates (x,y)
        if self.isidentity():
            return img
        if self.transfotype==transformationType.translation:
            # shift: takes transformation vector for coordinates (y,x)
            return scipy.ndimage.interpolation.shift(img,-self.cof[1::-1],cval = self.cval,**self.interpolationargs)
        else:
            if self.islinearidentity():
                # shift: takes transformation vector for coordinates (y,x)
                return scipy.ndimage.interpolation.shift(img,-self.cof[1::-1,2],cval = self.cval,**self.interpolationargs)
            elif self.isprojidentity():
                # affine_transform: takes change-of-frame matrix for coordinates (y,x)
                return scipy.ndimage.interpolation.affine_transform(img,self.cof[0:2,0:2].T,offset=self.cof[1::-1,2],cval = self.cval,**self.interpolationargs)
            else:
                raise NotImplementedError()

    def __repr__(self):
        if self.transfotype==transformationType.translation:
            return "tx = {}, ty = {}".format(self.cof[0],self.cof[1])
        else:
            return "tx = {}, ty = {}\nR={}\npx = {}, py = {}".format(self.cof[0,2],self.cof[1,2],self.cof[0:2,0:2],self.cof[2,0],self.cof[2,1])

    def __str__(self):
        return self.__repr__()

