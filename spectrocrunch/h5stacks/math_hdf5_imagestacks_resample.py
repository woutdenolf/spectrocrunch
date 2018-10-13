# -*- coding: utf-8 -*-
#
#   Copyright (C) 2017 European Synchrotron Radiation Facility, Grenoble, France
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
from scipy import interpolate

from ..io import nexus
from . import math_hdf5_imagestacks_copy as m_copy

class samplercls():

    def __init__(self,cval=np.nan,kind="linear"):
        self.newvalues = [None]*3 # 3x vector
        self.oldvalues = [None]*3 # 3x 3D dataset
        self.cval = cval
        self.kind = kind

    def getfirstaxis(self,func):
        for i in range(3):
            if func(self.newvalues[i]):
                return i
        return -1

    def getndims(self):
        return sum([v is not None for v in self.newvalues])

    def setencoderdim(self,axis,newvalues,oldvalues):
        self.newvalues[axis] = newvalues
        self.oldvalues[axis] = oldvalues

    def interpolate(self,datain,dataout):
        ndims = self.getndims()

        if ndims==0:
            dataout[:] = datain
        elif ndims==1:
            self._interpolate1D(datain,dataout,self.getfirstaxis(lambda v: v is not None))
        elif ndims==2:
            self._interpolate2D(datain,dataout,self.getfirstaxis(lambda v: v is None))
        else:
            self._interpolate3D(datain,dataout)

    def _interpolate1D(self,datain,dataout,axis):
        n0,n1,n2 = datain.shape

        xold = self.oldvalues[axis]
        xnew = self.newvalues[axis]

        if axis==0:
            for j in range(n1):
                for k in range(n2):
                    dataout[:,j,k] = self._interpolate1Ddo(xold[:,j,k], datain[:,j,k], xnew)
        elif axis==1:
            for j in range(n0):
                for k in range(n2):
                    dataout[j,:,k] = self._interpolate1Ddo(xold[j,:,k], datain[j,:,k], xnew)
        else:
            for j in range(n0):
                for k in range(n1):
                    dataout[j,k,:] = self._interpolate1Ddo(xold[j,k,:], datain[j,k,:], xnew)

    def _interpolate1Ddo(self,x,y,xnew):
        s = len(y)

        # Keep only data
        ind = np.isfinite(y)
        y = y[ind]

        if len(y)<2:
            return np.full(s,np.nan)

        x = x[ind]

        f = interpolate.interp1d(x,y,kind=self.kind, fill_value=self.cval, bounds_error=False, assume_sorted=False)
        return f(xnew)

    def _interpolate2D(self,datain,dataout,stackaxis):
        n0,n1,n2 = datain.shape

        if stackaxis==0:
            xold = self.oldvalues[2]
            yold = self.oldvalues[1]
            xnew = self.newvalues[2]
            ynew = self.newvalues[1]
            for i in range(n0):
                dataout[i,...] = self._interpolate2Ddo(xold[i,...],yold[i,...],datain[i,...],xnew,ynew)
        elif stackaxis==1:
            xold = self.oldvalues[2]
            yold = self.oldvalues[0]
            xnew = self.newvalues[2]
            ynew = self.newvalues[0]
            for i in range(n1):
                dataout[:,i,:] = self._interpolate2Ddo(xold[:,i,:],yold[:,i,:],datain[:,i,:],xnew,ynew)
        else:
            xold = self.oldvalues[1]
            yold = self.oldvalues[0]
            xnew = self.newvalues[1]
            ynew = self.newvalues[0]
            for i in range(n2):
                dataout[...,i] = self._interpolate2Ddo(xold[...,i],yold[...,i],datain[...,i],xnew,ynew)

    def _interpolate2Ddo(self,x,y,z,xnew,ynew):
        # f = interpolate.interp2d(x,y,z, kind=self.kind, fill_value=self.cval)
        # return f(xnew,ynew)

        s = z.shape

        # Keep only data
        z = z.ravel()
        ind = np.isfinite(z)
        z = z[ind]

        # No data left:
        if len(ind)<4:
            return np.full(s,np.nan)

        # Interpolate
        x = x.ravel()
        x = x[ind]
        y = y.ravel()
        y = y[ind]
        return interpolate.griddata((x,y),z, (xnew[None,:], ynew[:,None]), method=self.kind, fill_value=self.cval)

    def _interpolate3D(self,datain,dataout):
        raise NotImplementedError()

def evaluate_sliced(sampler,fin,varargs,retstacks):
    """ Interpolate data stacks
    Args:
        sampler(dict)
        fin(nexus.File)
        varargs(dict)
        fixedargs(dict)
        retstacks(list(h5py.Group))
    """
    for i in range(len(varargs)):
        grp = fin[varargs[i]]
        data = grp[grp.attrs["signal"]]
        dset = nexus.createNXdataSignal(retstacks[i],data=data,chunks = True)
        sampler.interpolate(data,dset)

def getencoderstacks(axes,encodermap,stacks):
    return [[s for s in stacks if encodermap[a["name"]]["encoder"] in s] if a["name"] in encodermap else [] for a in axes]

def getencoderresolution(x,y,axis):
    # res*x + offset = y
    #
    # (nX x 2) . (2 x K) = (nX x K)
    #
    #  x1 1   res   y1
    #  x2 1 . off = y2
    #  .  .         .
    #  xn 1         yn

    A = np.vstack([x, np.ones(len(x))]).T

    n0,n1,n2 = y.shape
    if axis==0:
        res = np.empty((n1,n2))
        off = np.empty((n1,n2))
        if n1<n2:
            for j in range(n1):
                res[j,:],off[j,:] =  np.linalg.lstsq(A,y[:,j,:])[0]
        else:
            for j in range(n2):
                res[:,j],off[:,j] =  np.linalg.lstsq(A,y[:,:,j])[0]

    elif axis==1:
        res = np.empty((n0,n2))
        off = np.empty((n0,n2))
        if n0<n2:
            for j in range(n0):
                res[j,:],off[j,:] =  np.linalg.lstsq(A,y[j,:,:])[0]
        else:
            for j in range(n2):
                res[:,j],off[:,j] =  np.linalg.lstsq(A,y[:,:,j].T)[0]
    else:
        res = np.empty((n0,n1))
        off = np.empty((n0,n1))
        if n0<n1:
            for j in range(n0):
                res[j,:],off[j,:] =  np.linalg.lstsq(A,y[j,:,:].T)[0]
        else:
            for j in range(n1):
                res[:,j],off[:,j] =  np.linalg.lstsq(A,y[:,j,:].T)[0]

    return np.median(res),np.median(off)

def calcoffset(diff,axis=None):
    return (np.max(diff,axis=axis)+np.min(diff,axis=axis))/2.

def getencoderoffset(xm,y,axis):
    
    n0,n1,n2 = y.shape
    if axis==0:
        off = np.empty((n1,n2))
        if n1<n2:
            for j in range(n1):
                off[j,:] = calcoffset(y[:,j,:]-xm[:,np.newaxis],axis=0)
        else:
            for j in range(n2):
                off[:,j] = calcoffset(y[:,:,j]-xm[:,np.newaxis],axis=0)
    elif axis==1:
        off = np.empty((n0,n2))
        if n0<n2:
            for j in range(n0):
                off[j,:] = calcoffset(y[j,:,:]-xm[:,np.newaxis],axis=0)
        else:
            for j in range(n2):
                off[:,j] = calcoffset(y[:,:,j]-xm[np.newaxis,:],axis=1)
    else:
        off = np.empty((n0,n1))
        if n0>n1:
            for j in range(n0):
                off[j,:] = calcoffset(y[j,:,:]-xm[np.newaxis,:],axis=1)
        else:
            for j in range(n1):
                off[:,j] = calcoffset(y[:,j,:]-xm[np.newaxis,:],axis=1)

    return np.median(off)

def analyseencoder(axes,encodermap,stacks,fin):

    # Select the encoder stack for each axis (can give no stack or more than one)
    encoderstackname = getencoderstacks(axes,encodermap,stacks)

    # Encoder = axis*res + offset
    sampler = samplercls()
    for axis in range(3):
        # Axis has encoder?
        if len(encoderstackname[axis])!=1:
            continue

        axisvalues = fin[axes[axis]["fullname"]]
        
        grp = fin[encoderstackname[axis][0]]
        encoderstack = grp[grp.attrs["signal"]]

        #res,off = getencoderresolution(axisvalues,encoderstack,axis)
        res = encodermap[axes[axis]["name"]]["resolution"]
        xm = axisvalues[:]*res
        off = getencoderoffset(xm,encoderstack,axis)

        sampler.setencoderdim(axis,xm+off,encoderstack)

    return sampler

def evaluate(operation,fin,varargs,retstacks,axes):
    """ Resample axes based on encoders
    Args:
        operation(dict)
        fin(nexus.File)
        varargs(dict)
        fixedargs(dict)
        retstacks(list(h5py.Group))
    """
    axesdata = []

    sampler = analyseencoder(axes,operation,varargs,fin)

    # Add (resampled) data to the NXdata group
    if sampler.getndims()==0:
        m_copy.evaluate(operation,fin,varargs,retstacks)
    else:
        evaluate_sliced(sampler,fin,varargs,retstacks)

        axesdata = [None]*len(axes)
        for i in range(len(axes)):
            axesdata[i] = fin[axes[i]["fullname"]][:]

    return axesdata

