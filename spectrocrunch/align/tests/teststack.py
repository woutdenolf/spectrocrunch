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
import scipy.ndimage

from ..types import transformationType

def makeGaussian(x,y,x0,y0,sx,sy,rho,A):
    num = (x-x0)**2/sx**2 - 2*rho/(sx*sy)*(x-x0)*(y-y0) + (y-y0)**2/sy**2
    denom = 2*(1-rho**2)
    return A*np.exp(-num/denom)#/(2*np.pi*sx*sy*np.sqrt(1-rho**2))

def makeGaussianAngle(x,y,x0,y0,sx,sy,angle,A):
    cosa = np.cos(angle)
    sina = np.sin(angle)
    sin2a = np.sin(2*angle)
    varx = sx**2
    vary = sy**2
    a = cosa**2/(2*varx)+sina**2/(2*vary)
    b = -sin2a/(4*varx)+sin2a/(4*vary)
    c = sina**2/(2*varx)+cosa**2/(2*vary)
    num = a*(x-x0)**2 - 2*b*(x-x0)*(y-y0) + c*(y-y0)**2
    return A*np.exp(-num)

def gettransformedimage(x,y,data,angle=False):
    ret = np.zeros(x.shape,dtype=float)
    if angle:
        proc = makeGaussianAngle
    else:
        proc = makeGaussian
    for x0,y0,sx,sy,rho,A in data:
        ret += proc(x,y,x0,y0,sx,sy,rho,A)
    #ret /= np.max(ret)
    return ret

#def getlargeimage(npeaks,nsigma,shape,subshape):
#    n = np.round(npeaks*min(shape)/min(subshape)).astype(np.int)
#    image = np.zeros(shape,dtype=np.float32)
#    np.random.seed(1)
#    x = shape[0]*np.random.random(n**2)
#    y = shape[1]*np.random.random(n**2)
#    image[x.astype(np.int), y.astype(np.int)] = 1
#    image = scipy.ndimage.gaussian_filter(image, sigma=min(subshape)/(npeaks*nsigma))
#    image /= np.max(image)

def translation(dx,dy):
    Mcoord = np.identity(3)
    Mcof = np.identity(3)
    Mcoord[0:2,2] = [dx,dy]
    Mcof[0:2,2] = [-dx,-dy]
    return Mcoord,Mcof

def rigid(a,tx,ty):
    if (a<0 or a>1):
        raise ValueError("A rigid transformation is expected.")
    Mcoord = np.identity(3)
    Mcof = np.identity(3)
    b = np.sqrt(1-a*a)
    Mcoord[0,0:3] = [a,-b,tx]
    Mcoord[1,0:3] = [b,a,ty]
    Mcof[0,0:3] = [a,b,-a*tx-b*ty]
    Mcof[1,0:3] = [-b,a,b*tx-a*ty]
    return Mcoord,Mcof

def similarity(a,b,tx,ty):
    Mcoord = np.identity(3)
    Mcof = np.identity(3)
    Mcoord[0,0:3] = [a,-b,tx]
    Mcoord[1,0:3] = [b,a,ty]
    s = np.float(a*a+b*b)
    Mcof[0,0:3] = [a/s,b/s,-(a*tx+b*ty)/s]
    Mcof[1,0:3] = [-b/s,a/s,(b*tx-a*ty)/s]
    return Mcoord,Mcof

def affine(a,b,c,d,tx,ty):
    Mcoord = np.identity(3)
    Mcof = np.identity(3)
    Mcoord[0,0:3] = [a,b,tx]
    Mcoord[1,0:3] = [c,d,ty]
    det = np.float(a*d-b*c)
    Mcof[0,0:3] = [d/det,-b/det,(d*tx-b*ty)/det]
    Mcof[1,0:3] = [-c/det,a/det,(c*tx-a*ty)/det]
    return Mcoord,Mcof

def homography(a,b,c,d,tx,ty,px,py):
    Mcoord = np.identity(3)
    Mcof = np.identity(3)
    Mcoord[0,0:3] = [a,b,tx]
    Mcoord[1,0:3] = [c,d,ty]
    Mcoord[2,0:3] = [px,py,1]
    det = np.float(a*d-a*py*ty-b*c+b*px*ty+c*py*tx-d*px*tx)
    Mcof[0,0:3] = [d-py*ty, py*tx-b, b*ty-d*tx]
    Mcof[1,0:3] = [px*ty-c, a-px*tx, c*tx-a*ty]
    Mcof[2,0:3] = [c*py-d*px, b*px-a*py,a*d-b*c]
    Mcof /= det
    return Mcoord,Mcof

def transformation(t,n):
    if t==transformationType.rigid:
        # total rotation is 40 degrees
        a = np.around(np.cos(40./180.*np.pi/(n-1)),3)
        #a = np.sqrt(3)/2 # between -1 and 1
    elif t==transformationType.similarity:
        a = np.around(np.cos(40./180.*np.pi/(n-1)),3)
        b = np.around(np.sin(40./180.*np.pi/(n-1)),3)
        a *= 1.1
        b *= 1.1
    else:
        a = np.around(np.cos(40./180.*np.pi/(n-1)),3)
        b = np.around(np.sin(40./180.*np.pi/(n-1)),3)
        a *= 1.1
        b *= 1.2
        c = -b*0.9
        d = a*0.9

    tx = 1.
    ty = -2.
    px = 0.001
    py = -0.001

    if t==transformationType.translation:
        return translation(tx,ty)
    elif t==transformationType.rigid:
        return rigid(a,tx,ty)
    elif t==transformationType.similarity:
        return similarity(a,b,tx,ty)
    elif t==transformationType.affine:
        return affine(a,b,c,d,tx,ty)
    elif t==transformationType.homography:
        return homography(a,b,c,d,tx,ty,px,py)
    else:
        raise NotImplementedError("This transformation to has not been implemented.")

def random(a,b,n):
    return a+(b-a)*np.random.random(n)

def teststack(transfotype,nimages = 5):
    """
    Returns:
        3-tuple: list of image stacks, change-of-coordinate matrix between the subsequent images, stack dimensions
    """

    # Transformation between images
    Mcoord,Mcof = transformation(transfotype,nimages)

    # Shape of a subimage
    subshape = (71,61)
    subxmin = -subshape[1]//2
    subymin = -subshape[0]//2
    subxmax = subshape[1]//2
    subymax = subshape[0]//2
    subxmax = subshape[1]-1
    subymax = subshape[0]-1
    subshape = (subymax-subymin+1,subxmax-subxmin+1)

    # Determine shape of large image (transform corners)
    xy = np.empty((3,4))
    xy[0,:] = [subxmin,subxmax,subxmin,subxmax]
    xy[1,:] = [subymin,subymin,subymax,subymax]
    xy[2,:] = [1,1,1,1]
    myminmax = np.append(np.min(xy,axis=1),np.max(xy,axis=1))

    for i in range(1,nimages):
        xy = np.dot(Mcoord,xy)
        xy[0,:] /= xy[2,:]
        xy[1,:] /= xy[2,:]
        myminmax[0:3] = np.minimum(myminmax[0:3],np.min(xy,axis=1))
        myminmax[3:] = np.maximum(myminmax[3:],np.max(xy,axis=1))
    xmin = int(myminmax[0])
    ymin = int(myminmax[1])
    xmax = int(myminmax[3])
    ymax = int(myminmax[4])

    # Gaussians in large image
    npeaks = 20
    sxy = min(subshape)/(1.5*npeaks)
    margin = int(10*sxy)
    xmin -= margin
    ymin -= margin
    xmax += margin
    ymax += margin
    shape = (ymax-ymin+1,xmax-xmin+1)

    n = int(np.round(npeaks*min(shape)/min(subshape)))
    np.random.seed(1)
    x0 = random(0,shape[1],n**2)+xmin
    y0 = random(0,shape[0],n**2)+ymin
    sx = random(sxy*0.8,sxy*1.2,n**2)
    sy = random(sxy*0.8,sxy*1.2,n**2)
    rho = random(0,0.2,n**2)
    A = random(1,2,n**2)
    data = zip(x0,y0,sx,sy,rho,A)

    # Plot full images
    if False:
        xv, yv = np.meshgrid(np.arange(xmin,xmax+1),np.arange(ymin,ymax+1))
        xv = xv.reshape((1,shape[0]*shape[1]))
        yv = yv.reshape((1,shape[0]*shape[1]))
        xy = np.vstack((xv,yv,np.ones_like(xv)))
        img = gettransformedimage(xv,yv,data).reshape(shape)
        import matplotlib.pyplot as plt
        plt.figure(2)
        plt.subplot(111)
        plt.imshow(img,origin='lower',interpolation='nearest')
        plt.show()

    # Stack of transformed subimages
    ret = np.empty(subshape+(nimages,),dtype=np.float32)
    xv, yv = np.meshgrid(np.arange(subxmin,subxmax+1),np.arange(subymin,subymax+1))
    xv = xv.reshape((1,subshape[0]*subshape[1]))
    yv = yv.reshape((1,subshape[0]*subshape[1]))
    xy = np.vstack((xv,yv,np.ones_like(xv)))
    ret[...,0] = gettransformedimage(xv,yv,data).reshape(subshape)
    for i in range(1,nimages):
        xy = np.dot(Mcof,xy) # coordinates from new frame to old frame
        xy[0,:] /= xy[2,:]
        xy[1,:] /= xy[2,:]
        ret[...,i] = gettransformedimage(xy[0,:],xy[1,:],data).reshape(subshape)

    # Relative change-of-frame in subimage pixel frame
    C = np.identity(3,dtype=Mcof.dtype)
    Cinv = np.identity(3,dtype=Mcof.dtype)
    C[0:2,2] = [subxmin,subymin] # image pixel frame to subimage pixel frame
    Cinv[0:2,2] = -C[0:2,2]
    Mcof = np.dot(np.dot(Cinv,Mcof),C)
    Mcoord = np.dot(np.dot(Cinv,Mcoord),C)

    return ([ret,ret],Mcoord,2) # Mcoord: change-of-frame matrix of the back-transformation

