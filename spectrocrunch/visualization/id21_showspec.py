# -*- coding: utf-8 -*-
#
#   Copyright (C) 2016 European Synchrotron Radiation Facility, Grenoble, France
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

from spectrocrunch.io.spec import spec
import h5py
import numpy as np
import fabio
import os
from matplotlib import cm
import matplotlib.pyplot as plt

class shape_object(object):
    
    def __init__(self,noborder=False,cmap="jet",lo=0.05,hi=0.85):
        self.noborder=noborder
        self.cmap = cmap
        self.lo = lo
        self.hi = hi
        self.img = None

    def unloadimage(self):
        self.img = None

    def setimage(self,images):
        images = [np.zeros(s) if img is None else img for img in images]
        self.img = np.dstack(images)
        if len(self.img)==2:
            self.img = self.img[...,np.newaxis]

    def get2ddims(self):
        dim1 = np.array(self.dim1)
        dim2 = np.array(self.dim2)
        shape = 2
        if len(dim1)==1:
            dim1 = np.array([dim1,dim1])
            shape -= 1
        if len(dim2)==1:
            dim2 = np.array([dim2,dim2])
            shape -= 1
        return dim1.tolist(),dim2.tolist(),shape

    def getimage(self,transform,origin):
        # Image
        img = self.img
        if img is None:
            nchan = 0
            n1 = 0
            n2 = 0
        else:
            n1,n2,nchan = self.img.shape
            # Transform
            for c in transform:
                if c=='t': # Transpose
                    img = img.transpose((1,0,2))
                elif c=='v': # Flip vertical
                    img = img[::-1,:,:]
                elif c=='h': # Flip horizontal
                    img = img[:,::-1,:]

        # Plot range
        dim1,dim2,shape = self.get2ddims()
        if n1<=1:
            d1 = 0
        else:
            d1 = (dim1[1]-dim1[0])/(2.*n1-2)
        if n2<=1:
            d2 = 0
        else:
            d2 = (dim2[1]-dim2[0])/(2.*n2-2)
        
        extent = (dim2[0]-d2-origin[1],dim2[1]+d2-origin[1],dim1[0]-d1-origin[0],dim1[1]+d1-origin[0])
   
        for c in transform:
            if c=='t': # Transpose
                extent = (extent[2],extent[3],extent[0],extent[1])
            elif c=='v': # Flip vertical
                extent = (extent[0],extent[1],extent[3],extent[2])
            elif c=='h': # Flip horizontal
                extent = (extent[1],extent[0],extent[2],extent[3])

        return img,extent,nchan,shape

    def plot(self,ax,origin,transform=""):
        # Load the image
        #self.loadimage()

        # Get image
        img,extent,nchan,shape = self.getimage(transform,origin)

        # Images scaling
        for i in range(nchan):
            mi = np.min(img[...,i])
            ma = np.max(img[...,i])
            mi += (ma-mi)*self.lo
            ma -= (ma-mi)*self.hi
            img[...,i] -= mi
            img[...,i] /= ma
            img[...,i] = np.clip(img[...,i],0,1)

        # Plot image (if any)
        if nchan==1:
            ax.imshow(np.squeeze(img),extent=extent,origin='lower',interpolation='nearest',aspect=1,cmap=cm.get_cmap(self.cmap))
        elif nchan==2:
            ax.imshow(np.concatenate((img,np.zeros(img.shape[0:2])[...,np.newaxis]),axis=2),extent=extent,origin='lower',interpolation='nearest',aspect=1)
        elif nchan==3:
            ax.imshow(img,extent=extent,origin='lower',interpolation='nearest',aspect=1)

        # Plot border (or line or point)
        if not self.noborder:
            if shape==0:
                ax.scatter(extent[0], extent[2], marker='o')
            elif shape==1:
                if extent[0]==extent[1]:
                    ax.plot([extent[0],extent[0]], [extent[2],extent[3]])
                else:
                    ax.plot([extent[0],extent[1]], [extent[2],extent[2]])
            else:
                ax.plot([extent[0],extent[0],extent[1],extent[1],extent[0]],\
                        [extent[2],extent[3],extent[3],extent[2],extent[2]])

        # Free memory image
        #self.unloadimage()

class shape_hdf5(shape_object):

    def __init__(self,filename,subpaths,subindex,**kwargs):
        shape_object.__init__(self,**kwargs)
        self.filename = filename
        self.subpaths = subpaths
        self.subindex = subindex
        self.loadimage()

    def loadimage(self):
        oh5 = h5py.File(self.filename)

        # Prepare global coordinates
        dim1off = 0.
        dim1name = "samz"
        dim1mult = 1
        dim2off = 0.
        dim2name = "samy"
        dim2mult = 1

        ocoord = oh5["coordinates"]
        for f in ocoord:
            if f == "samz":
                v = np.atleast_1d([ocoord[f].value])
                if len(v)==1:
                    dim1off = v[0]*1000
                else:
                    dim1off = v[self.subindex]*1000

                dim1name = "sampz"
                dim1mult = 1
            if f == "sampz":
                v = np.atleast_1d([ocoord[f].value])
                if len(v)==1:
                    dim1off = v[0]
                else:
                    dim1off = v[self.subindex]

                dim1name = "samz"
                dim1mult = 1000
            if f == "samy":
                v = np.atleast_1d([ocoord[f].value])
                print(v)
                if len(v)==1:
                    dim2off = v[0]*1000
                else:
                    dim2off = v[self.subindex]*1000

                dim2name = "sampy"
                dim2mult = 1
            if f == "sampy":
                v = np.atleast_1d([ocoord[f].value])
                if len(v)==1:
                    dim2off = v[0]
                else:
                    dim2off = v[self.subindex]

                dim2name = "samy"
                dim2mult = 1000

        # Get image with axes in micron
        n = len(self.subpaths)
        images = [None]*n
        for i in range(n):
            ogrp = oh5[self.subpaths[i]]
            odset = ogrp[ogrp.attrs["signal"]]
            self.dim1 = dim1off + ogrp[dim1name].value[[0,-1]]*dim1mult
            self.dim2 = dim2off + ogrp[dim2name].value[[0,-1]]*dim2mult

            tmp = ogrp.attrs["axes"].split(":")
            idim1 = tmp.index(dim1name)
            idim2 = tmp.index(dim2name)
            if idim2!=0 and idim1!=0:
                img = odset[self.subindex,...]
            elif idim2!=1 and idim1!=1:
                img = odset[:,self.subindex,:]
            else:
                img = odset[...,self.subindex]
            if idim1 > idim2:
                img = img.T
            images[i] = [img]
        self.setimage(images)

        oh5.close()

class shape_spec(shape_object):
    
    def __init__(self,specfile,scannumber,labels=[],olddir="",newdir="",**kwargs):
        shape_object.__init__(self,**kwargs)
        self.specfile = specfile
        self.scannumber = scannumber
        self.labels = labels
        self.olddir = olddir
        self.newdir = newdir
        self.loadimage()

    def loadimage(self):
        motors = ["samz","sampz","samy","sampy","Energy MONO"]

        f = spec(self.specfile)

        p = f.getdimensions(self.scannumber,motors)
        self.dim1 = p["samz"]*1000 + p["sampz"]
        self.dim2 = p["samy"]*1000 + p["sampy"]

        print("Energy = {} keV".format(p["Energy MONO"]))
        
        # Get images
        n = len(self.labels)
        if n!=0:
            xia = f.getxialocation(self.scannumber)
            xia["DIRECTORY"] = xia["DIRECTORY"].replace(self.olddir,self.newdir)
            xia["ZAP SCAN NUMBER"] = int(xia["ZAP SCAN NUMBER"])
            xia["ZAP IMAGE NUMBER"] = int(xia["ZAP IMAGE NUMBER"])

            images = [None]*n
            s = None
            for i in range(n):
                filename = "%s/%s_%s_%04d_%04d.edf"%(\
                            xia["DIRECTORY"],xia["RADIX"],\
                            self.labels[i],xia["ZAP SCAN NUMBER"],\
                            xia["ZAP IMAGE NUMBER"])
                if os.path.isfile(os.path.join(filename)):
                    images[i] = fabio.open(filename).data
                    s = images[i].shape
                    print("  {}".format(filename))
                else:
                    print("  File not found: {}".format(filename))
            if s is not None:
                self.setimage(images)
            
def totaldims(assocs):
    dim1 = []
    dim2 = []
    for a in assocs:
        tmp1,tmp2,_ = a.get2ddims()
        dim1 += tmp1
        dim2 += tmp2
    return [min(dim1),max(dim1),min(dim2),max(dim2)]

def plot(lst):
    f,ax = plt.subplots(1)

    dim1min,dim1max,dim2min,dim2max = totaldims(lst)
    origin = [dim1min,dim2min]
    xlimits = [0,abs(dim2max-dim2min)]
    ylimits = [0,abs(dim1max-dim1min)]

    for l in lst:
        l.plot(ax,origin)

    ax.set_xlabel('X ($\mu$m)')
    ax.set_ylabel('Y ($\mu$m)')

    ax.set_xlim(xlimits[0],xlimits[1])
    ax.set_ylim(ylimits[0],ylimits[1])

    plt.show()





