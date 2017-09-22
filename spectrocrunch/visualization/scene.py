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

import matplotlib.pyplot as plt
import matplotlib.colors as pltcolors
import numpy as np

from ..common import listtools
from ..common import instance
from ..common.hashable import Hashable
from collections import OrderedDict
from .. import ureg

def ColorNorm(name,*args):
    if name=="LogNorm":
        norm = lambda **kwargs: pltcolors.LogNorm(**kwargs)
    elif name=="PowerNorm":
        norm = lambda **kwargs: pltcolors.PowerNorm(*args,**kwargs)
    elif name=="SymLogNorm":
        norm = lambda **kwargs: pltcolors.SymLogNorm(*args,**kwargs)
    else:
        norm = pltcolors.Normalize
    return norm
    
class Object(object):
    
    def dimindex(self,dim):
        if self.transposed:
            return 1-dim
        else:
            return dim

    @property
    def xindex(self):
        return self.dimindex(1)

    @property
    def yindex(self):
        return self.dimindex(0)

class Scene(Object):
    """Each scene can have a number of registered items with associalted plot settings.
    """

    def __init__(self,ax=None,unit0=ureg.dimensionless,unit1=ureg.dimensionless,title=""):
        self.ax = ax
        self._items = OrderedDict()
        self.transposed = False
        self.axlabels = ["Dim1","Dim0"]
        self.cmap = plt.get_cmap('jet')
        self.decreasing = [False,False]
        self.title = title
        self.aspectratio = None
        self.units = [unit0,unit1]
        self.dataoffset = [ureg.Quantity(0,unit0),ureg.Quantity(0,unit1)]
        self.datascale = [1,1]
        
    def dimquantity(self,v,dim):
        if isinstance(v,ureg.Quantity):
            return v
        else:
            return ureg.Quantity(v,self.units[dim])
    
    def q0(self,v):
        return self.dimquantity(v,0)
    
    def q1(self,v):
        return self.dimquantity(v,1)
        
    def magnitude(self,v,dim):
        unit = self.units[dim]
        if instance.isarray(v):
            return type(v)(a.to(unit).magnitude for a in v)
        else:
            return v.to(unit).magnitude
    
    @property
    def xunits(self):
        return self.units[self.xindex]
    
    @property
    def yunits(self):
        return self.units[self.yindex]
     
    def xmagnitude(self,v):
        unit = self.xunits
        if instance.isarray(v):
            return type(v)(a.to(unit).magnitude for a in v)
        else:
            return v.to(unit).magnitude
            
    def ymagnitude(self,v):
        unit = self.yunits
        if instance.isarray(v):
            return type(v)(a.to(unit).magnitude for a in v)
        else:
            return v.to(unit).magnitude
    
    @property
    def xlabel(self):
        xunit = self.xunits
        if xunit == ureg.dimensionless:
            return self.axlabels[self.xindex]
        else:
            return "{} ({})".format(self.axlabels[self.xindex],xunit)
    
    @property
    def ylabel(self):
        yunit = self.yunits
        if yunit == ureg.dimensionless:
            return self.axlabels[self.yindex]
        else:
            return "{} ({})".format(self.axlabels[self.yindex],yunit)
        
    def setaxes(self,ax):
        self.ax = ax
        
    def transpose(self,tr):
        self.transposed = tr

    def register(self,item):
        item.addscene(self)
        self._items[item] = item.defaultsettings()

    def getitemsettings(self,item):
        if item in self._items:
            return self._items[item]
        else:
            return {}

    def resetdatarange(self):
        self.dataoffset = [self.q0(0),self.q1(1)]
        self.datascale = [1,1]
    
    @property
    def items(self):
        for item in self._items:
            item.selectscene(self)
        return self._items
        
    def datarange(self,dim):
        ran = [r for item in self.items for r in item.datarange(dim)]
        return min(ran),max(ran)

    def setdatarange(self,dim,ran):
        self.dataoffset[dim] = self.dimquantity(0,dim)
        self.datascale[dim] = 1
        ranorg = sorted(self.datarange(dim))
        scale = (ran[1]-ran[0])/(ranorg[1]-ranorg[0])
        off = scale*ranorg[0]-ran[0]
        self.dataoffset[dim] = off
        self.datascale[dim] = scale
    
    def origin(self,dim,off):
        self.dataoffset[dim] = off-self.dataoffset[dim]
    
    def originzero(self):
        self.origin(0,self.q0(0))
        self.origin(1,self.q1(0))
    
    def datatransform(self,arr,dim):
        if instance.isarray(arr):
            return type(arr)(x*self.datascale[dim] - self.dataoffset[dim] for x in arr)
        else:
            return arr*self.datascale[dim] - self.dataoffset[dim]
            
    @property
    def xlim(self):
        xlim = [x for item in self.items for x in item.xlim]
        if self.decreasing[self.xindex]:
            return max(xlim),min(xlim)
        else:
            return min(xlim),max(xlim)

    @property
    def ylim(self):
        ylim = [y for item in self.items for y in item.ylim]
        if self.decreasing[self.yindex]:
            return max(ylim),min(ylim)
        else:
            return min(ylim),max(ylim)

    @property
    def aspect(self):
        if self.aspectratio is None:
            return None
        ylim = self.ylim
        xlim = self.xlim
        return abs(xlim[1]-xlim[0])/abs(ylim[1]-ylim[0])*self.aspectratio
        
    def crop(self):
        lim = self.xmagnitude(self.xlim)
        self.ax.set_xlim(lim[0],lim[1])
        lim = self.ymagnitude(self.ylim)
        self.ax.set_ylim(lim[0],lim[1])

    def update(self,**kwargs):
        for item in self.items:
            item.update(**kwargs)
        self.crop()
        
        self.ax.set_title(self.title)
        self.ax.set_xlabel(self.xlabel)
        self.ax.set_ylabel(self.ylabel)
    
    @property
    def vminmax(self):
        vmin = None
        vmax = None
        
        single = lambda v: (np.nanmin(v[0]),np.nanmax(v[1]))
        
        tmp = [single(item.vminmax) for item in self.items if hasattr(item,'vminmax')]
        
        if len(tmp)==0:
            return None,None
        else:
            vmin,vmax = zip(*tmp)
            return min(vmin),max(vmax)
    
    def scale(self,vmin=None,vmax=None):
        self.set_setting("vmin",vmin)
        self.set_setting("vmax",vmax)
        
    def set_setting(self,attr,value):
        for item,settings in self._items.items():
            if attr in settings:
                settings[attr] = value
            
class Item(Hashable,Object):
    """Each item can be registered with multiple scenes. Data is owned by the Item instance,
       plot settings are owned by the scenes.
    """
    
    def __init__(self,dim0name="Dim0",dim1name="Dim1"):
        self._scenes = []
        self._plotobjs = []
        self._sceneindex = -1
        self.dim0name = dim0name
        self.dim1name = dim1name
        
    def useaxesnames(self):
        self.scene.axlabels = [self.dim0name,self.dim1name]
    
    def defaultsettings(self):
        return {}
    
    def set_setting(self,attr,value):
        settings = self.scene.getitemsettings(self)
        if attr in settings:
            settings[attr] = value
    
    def get_setting(self, key):
        settings = self.scene.getitemsettings(self)
        return settings[key]
    
    @property
    def scene(self):
        if len(self._scenes)==0:
            raise RuntimeError("This object is not registered with any scene")
        return self._scenes[self._sceneindex]
    
    def addscene(self,s):
        self._scenes.append(s)
        
    def selectscene(self,s):
        try:
            self._sceneindex = self._scenes.index(s)
        except:
            raise RuntimeError("This object is not registered with scene {}".format(s))
    
    @property
    def transposed(self):
        return self.scene.transposed
       
    @property
    def plotobjs(self):
        for ind,o in enumerate(self._plotobjs):
            for i in o:
                if i is not None:
                    if self.scene.ax == i.axes:
                        return o
        return []
    
    def addobjs(self,o):
        if isinstance(o,list):
            self._plotobjs.append(o)
        else:
            self._plotobjs.append([o])
    
    def datarange(self,dim):
        raise NotImplementedError("Item is an abstract class")
    
    @property
    def xlim(self):
        return self.datarange(self.xindex)
    
    @property
    def ylim(self):
        return self.datarange(self.yindex)
    
    def _stringrepr(self):
        return "{} {}".format(type(self).__name__,id(self))
    
class Image(Item):

    def __init__(self,img,lim0=None,lim1=None,**kwargs):
        self.img = img
        
        if lim0 is None:
            lim0 = [0,img.shape[0]-1]
        else:
            lim0 = [lim0[0],lim0[-1]]
            
        if lim1 is None:
            lim1 = [0,img.shape[1]-1]
        else:  
            lim1 = [lim1[0],lim1[-1]]
    
        self.lim = [lim0,lim1]
        
        super(Image,self).__init__(**kwargs)
    
    @property
    def vminmax(self):
        if self.rgb:
            n = min(3,self.img.ndim)
            vmin = np.asarray([np.nanmin(self.img[...,i]) for i in range(n)])
            vmax = np.asarray([np.nanmax(self.img[...,i]) for i in range(n)])
        else:
            vmin = np.nanmin(self.img)
            vmax = np.nanmax(self.img)
        return vmin,vmax
    
    def scale(self,vmin=None,vmax=None):
        self.set_setting("vmin",vmin)
        self.set_setting("vmax",vmax)
    
    def selfscale(self,mvmin=0,mvmax=1):
        vmin,vmax = self.vminmax
        d = vmax-vmin
        if instance.isarray(d):
            if instance.isarray(mvmin):
                mvmin = np.asarray(mvmin)
            if instance.isarray(mvmax):
                mvmax = np.asarray(mvmax)
        vmax = vmin + d*mvmax
        vmin = vmin + d*mvmin
        self.scale(vmin=vmin,vmax=vmax)
    
    @staticmethod
    def defaultsettings():
        return {"cmap":None,\
                  "aspect":None,\
                  "plotborder":True,\
                  "color":None,\
                  "plotimage":True,\
                  "linewidth":2,\
                  "alpha":1,\
                  "vmin":None,\
                  "vmax":None,\
                  "cnorm":None}
        
    def datarange(self,dim,border=False):
        lim = self.scene.datatransform(self.lim[dim],dim)
        if border:
            d = (lim[1]-lim[0])/(2.*(self.img.shape[dim]-1))
            return lim[0]-d,lim[1]+d
        else:
            return lim
    
    @property
    def xlim(self):
        return self.datarange(self.xindex,border=True)
    
    @property
    def ylim(self):
        return self.datarange(self.yindex,border=True)

    def bordercoord(self):
        x = sorted(self.xlim)
        y = sorted(self.ylim)
        x = [x[0],x[1],x[1],x[0],x[0]]
        y = [y[0],y[0],y[1],y[1],y[0]]
        return x,y
        
    @property
    def image(self):
        if self.transposed:
            return np.swapaxes(self.img,0,1)
        else:
            return self.img

    @property
    def rgb(self):
        return self.img.ndim>2
        
    def update(self):  
        scene = self.scene
        
        settings = scene.getitemsettings(self)
        if settings["cmap"] is None:
            cmap = scene.cmap
        else:
            cmap = settings["cmap"]
        if settings["aspect"] is None:
            aspect = scene.aspect
        else:
            aspect = settings["aspect"]
        alpha = settings["alpha"]
        
        vmin,vmax = scene.vminmax
        if settings["vmin"] is not None:
            vmin = settings["vmin"]   
        if settings["vmax"] is not None:
            vmax = settings["vmax"]
        image = self.image
            
        if self.rgb:
            if not instance.isarray(vmin):
                vmin = [vmin]*self.img.ndim
            if not instance.isarray(vmax):
                vmax = [vmax]*self.img.ndim

            print ""
            print vmin,vmax

            for i in range(self.img.ndim):
                if settings["cnorm"] is None:
                    norm = pltcolors.Normalize(vmin=vmin[i],vmax=vmax[i],clip=True)
                else:
                    norm = settings["cnorm"](vmin=vmin[i],vmax=vmax[i],clip=True)
                    
                print np.nanmin(image[...,i]),np.nanmax(image[...,i])
                image[...,i] = norm(image[...,i]).data
                #print np.nanmin(imgi),np.nanmax(imgi)
                
                #if vmin[i]==vmax[i]:
                #    image[...,i] = imgi
                #else:
                #    image[...,i] = (imgi-vmin[i])/(vmax[i]-vmin[i])
                
                print np.nanmin(image[...,i]),np.nanmax(image[...,i])
                
            norm = None
        else:
            if settings["cnorm"] is None:
                norm = pltcolors.Normalize(vmin=vmin,vmax=vmax)
            else:
                norm = settings["cnorm"](vmin=vmin,vmax=vmax)
            
        
        
        extent = scene.xmagnitude(self.xlim)+scene.ymagnitude(self.ylim)
        
        x,y = self.bordercoord()
        x = scene.xmagnitude(x)
        y = scene.ymagnitude(y)
        
        o = self.plotobjs
        if len(o)==0:
            # matplotlib.image.AxesImage
            o1 = scene.ax.imshow(image,\
                        interpolation = 'nearest',\
                        extent = extent,\
                        cmap = cmap,\
                        origin = 'lower',
                        norm = norm,\
                        alpha = alpha,\
                        aspect = aspect)
            o1.set_visible(settings["plotimage"])
            
            o2 = scene.ax.plot(x, y, color=settings["color"], linewidth=settings["linewidth"])[0]
            o2.set_visible(settings["plotborder"])
            settings["color"] = o2.get_color()
            
            self.addobjs([o1,o2])
        else:
            o[0].set_data(image)
            plt.setp(o[0], interpolation = 'nearest',\
                        extent = extent,\
                        norm = norm,\
                        alpha = alpha,\
                        cmap = cmap)
            o[0].set_visible(settings["plotimage"])
            
            o[1].set_data(x,y)
            plt.setp(o[1], color=settings["color"], linewidth=settings["linewidth"])
            o[1].set_visible(settings["plotborder"])
            settings["color"] = o[1].get_color()
            
class Polygon(Item):

    def __init__(self,pdim0,pdim1,scatter=False,closed=True,**kwargs):
        self._pdim = [pdim0,pdim1]
        self.scatter = scatter
        self.closed = closed

        super(Polygon,self).__init__(**kwargs)
    
    @staticmethod
    def defaultsettings():
        return {"color":None,"linewidth":2}
    
    def datarange(self,dim):
        ran = self.scene.datatransform(self._pdim[dim],dim)
        return min(ran),max(ran)
    
    def pdim(self,dim):
        idim = self.dimindex(dim)
        return self.scene.datatransform(self._pdim[idim],idim)
    
    def update(self):
        scene = self.scene
        settings = scene.getitemsettings(self)
        
        x = scene.xmagnitude(self.pdim(1))
        y = scene.ymagnitude(self.pdim(0))
        
        bscatter = not instance.isarray(x) or self.scatter
        if not bscatter and self.closed:
            x.append(x[0])
            y.append(y[0])
        
        o = self.plotobjs
        if len(o)==0:
            # matplotlib.lines.Line2D
            if bscatter:
                o = [scene.ax.plot(x, y, marker='o', linestyle = 'None', color=settings["color"], linewidth=settings["linewidth"])[0]]
            else:
                o = [scene.ax.plot(x, y, color=settings["color"], linewidth=settings["linewidth"])[0]]
            self.addobjs(o)
        else:
            o[0].set_data(x,y)
            plt.setp(o[0], color=settings["color"], linewidth=settings["linewidth"])
            
        settings["color"] = o[0].get_color()
        
