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

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as pltcolors

import numpy as np
import logging

from ..common import listtools
from ..common import instance
from ..common.hashable import Hashable
from collections import OrderedDict
from .. import ureg

logger = logging.getLogger(__name__)

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
    
    def showaxes(self):
        self.ax.set_axis_on()
        self.ax.get_xaxis().set_visible(True)
        self.ax.get_yaxis().set_visible(True)
    
    def hideaxes(self):
        self.ax.set_axis_off()
        self.ax.get_xaxis().set_visible(False)
        self.ax.get_yaxis().set_visible(False)
     
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
    
    def flipdim(self,dim):
        self.setdatarange(dim,self.datarange(dim)[::-1])

    def origin(self,dim,off):
        self.dataoffset[dim] = off-self.dataoffset[dim]
    
    def originreset(self):
        self.origin(0,self.q0(0))
        self.origin(1,self.q1(0))
    
    def originzero(self):
        self.origin(0,self.datarange(0)[0])
        self.origin(1,self.datarange(1)[0])
        
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
        dx = abs(xlim[1]-xlim[0])
        dy = abs(ylim[1]-ylim[0])
        return dx/dy*self.aspectratio
        
    def crop(self):
        lim = self.xmagnitude(self.xlim)
        self.ax.set_xlim(lim[0],lim[1])
        #self.ax.set_xbound(lim[0],lim[1])
        lim = self.ymagnitude(self.ylim)
        self.ax.set_ylim(lim[0],lim[1])
        #self.ax.set_ybound(lim[0],lim[1])
        
    def update(self,**kwargs):
        for item in self.items:
            item.update(**kwargs)
        self.crop()
        
        if self.title is not None:
            self.ax.set_title(self.title)
        if self.xlabel is not None:
            self.ax.set_xlabel(self.xlabel)
        if self.ylabel is not None:
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
    
    def __init__(self,dim0name="Dim0",dim1name="Dim1",scene=None,**kwargs):
        self._scenes = []
        self._plotobjs = [] # list of lists, each list belonging to a scene
        self._sceneindex = -1
        self.dim0name = dim0name
        self.dim1name = dim1name
        
        if scene is not None:
            self.register(scene)
        
        try:
            self.set_settings(kwargs)
        except RuntimeError:
            logger.warning("Item settings are not applied (provide a scene)")
    
    def register(self,scene):
        scene.register(self)
    
    def useaxesnames(self):
        self.scene.axlabels = [self.dim0name,self.dim1name]
    
    def defaultsettings(self):
        return {}
    
    def set_settings(self,settings):
        settings = self.scene.getitemsettings(self)
        for k,v in settings.items():
            settings[k] = v

    def get_settings(self, keys):
        settings = self.scene.getitemsettings(self)
        return {k:settings[k] for k in keys}
        
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
        """Objects belonging to the active scene
        """
        for objs in self._plotobjs:
            for o in objs:
                if o is not None:
                    if self.scene.ax == o.axes:
                        return objs
        return []
    
    def addobjs(self,objs):
        """Objects belonging to the active scene
        """
        self.delobjs()
        if isinstance(objs,list):
            self._plotobjs.append(objs)
        else:
            self._plotobjs.append([objs])
    
    def delobjs(self):
        objs = self.plotobjs
        if len(objs)>0:
            for o in objs:
                if o is not None:
                    o.remove()
            self._plotobjs.remove(objs)
            

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
        self._img = img
        
        if lim0 is None:
            lim0 = [0,img.shape[0]-1]
        else:
            lim0 = [lim0[0],lim0[-1]]
            
        if lim1 is None:
            lim1 = [0,img.shape[1]-1]
        else:  
            lim1 = [lim1[0],lim1[-1]]
    
        self.lim = [lim0,lim1]

        self.cax = None
        
        super(Image,self).__init__(**kwargs)
    
    @property
    def image(self):
        if self.transposed:
            return np.swapaxes(self.img,0,1)
        else:
            return self.img
            
    @property
    def img(self):
        if self.rgb:
            img = np.zeros(self.shape,self.dtype)
            for i,j in enumerate(self.channels):
                img[...,i] = self._img[...,j]
            return img
        else:
            return self._img.copy()

    @property
    def rgb(self):
        return self._img.ndim==3
        
    @property
    def nimg(self):
        if self.rgb:
            return min(self._img.shape[2],3)
        else:
            return 1

    @property
    def shape(self):
        if self.rgb:
            return self._img.shape[0],self._img.shape[1],3
        else:
            return self._img.shape[0],self._img.shape[1]
    
    @property
    def dtype(self):
        return self._img.dtype
        
    @property
    def channels(self):
        settings = self.scene.getitemsettings(self)
        ind = settings["channels"]
        if ind is None:
            return range(self.nimg)
        if len(ind)>3:
            return ind[0:3]
        return ind

    @property
    def vminmax(self):
        img = self.img
        if self.rgb:
            vmin = np.asarray([np.nanmin(img[...,i]) for i in range(self.nimg)])
            vmax = np.asarray([np.nanmax(img[...,i]) for i in range(self.nimg)])
        else:
            vmin = np.nanmin(img)
            vmax = np.nanmax(img)
        return vmin,vmax
    
    def scale(self,vmin=None,vmax=None):
        self.set_setting("vmin",vmin)
        self.set_setting("vmax",vmax)
    
    def selfscale(self,pmin=0,pmax=1):
        vmin,vmax = self.vminmax
        d = vmax-vmin
        if instance.isarray(d):
            if instance.isarray(pmin):
                pmin = np.asarray(pmin)
            if instance.isarray(pmax):
                pmax = np.asarray(pmax)
        vmax = vmin + d*pmax
        vmin = vmin + d*pmin
        self.scale(vmin=vmin,vmax=vmax)
    
    @staticmethod
    def defaultsettings():
        return {"cmap":None,\
                  "aspect":None,\
                  "border":False,\
                  "color":None,\
                  "image":True,\
                  "linewidth":2,\
                  "alpha":1,\
                  "vmin":None,\
                  "vmax":None,\
                  "cnorm":None,\
                  "channels":None}
        
    def datarange(self,dim,border=False):
        lim = self.scene.datatransform(self.lim[dim],dim)
        if border:
            d = (lim[1]-lim[0])/(2.*(self.shape[dim]-1))
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
    def aspect(self):
        scene = self.scene
        settings = scene.getitemsettings(self)
        if settings["aspect"] is None:
            return scene.aspect
        else:
            return settings["aspect"]
            
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
                vmin = [vmin]*self.nimg
            if not instance.isarray(vmax):
                vmax = [vmax]*self.nimg
            for i in range(self.nimg):
                if settings["cnorm"] is None:
                    norm = pltcolors.Normalize(vmin=vmin[i],vmax=vmax[i],clip=True)
                else:
                    norm = settings["cnorm"](vmin=vmin[i],vmax=vmax[i],clip=True)
                    
                image[...,i] = norm(image[...,i]).data
                
                #print np.nanmin(imgi),np.nanmax(imgi)
                
                #if vmin[i]==vmax[i]:
                #    image[...,i] = imgi
                #else:
                #    image[...,i] = (imgi-vmin[i])/(vmax[i]-vmin[i])
                
                #print np.nanmin(image[...,i]),np.nanmax(image[...,i])
                
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
            self.cax = o1 = scene.ax.imshow(image,\
                        interpolation = 'nearest',\
                        extent = extent,\
                        cmap = cmap,\
                        origin = 'lower',
                        norm = norm,\
                        alpha = alpha,\
                        aspect = aspect)
            o1.set_visible(settings["image"])
            
            o2 = scene.ax.plot(x, y, color=settings["color"], linewidth=settings["linewidth"])[0]
            o2.set_visible(settings["border"])
            settings["color"] = o2.get_color()
            
            self.addobjs([o1,o2])
        else:
            o[0].set_data(image)
            plt.setp(o[0], interpolation = 'nearest',\
                        extent = extent,\
                        norm = norm,\
                        alpha = alpha,\
                        cmap = cmap)
            o[0].set_visible(settings["image"])
            
            o[1].set_data(x,y)
            plt.setp(o[1], color=settings["color"], linewidth=settings["linewidth"])
            o[1].set_visible(settings["border"])
            settings["color"] = o[1].get_color()
            
class Polyline(Item):

    def __init__(self,pdim0,pdim1,**kwargs):
        self._pdim = [pdim0,pdim1]
        super(Polyline,self).__init__(**kwargs)
    
    @staticmethod
    def defaultsettings():
        return {"color":None,"linewidth":2,"scatter":False,"closed":True}
    
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
        
        bscatter = not instance.isarray(x) or settings["scatter"]
        if not bscatter and settings["closed"]:
            x = list(x)
            y = list(y)
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

class Polygon(Item):

    def __init__(self,pdim0,pdim1,**kwargs):
        self._pdim = [pdim0,pdim1]

        super(Polygon,self).__init__(**kwargs)
    
    @staticmethod
    def defaultsettings():
        return {"color":None,"linewidth":2,"closed":True}
    
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
        
        if settings["closed"]:
            x = list(x)
            y = list(y)
            x.append(x[0])
            y.append(y[0])
        
        o = [matplotlib.patches.Polygon(zip(x, y), color=settings["color"], linewidth=settings["linewidth"])]
        scene.ax.add_patch(o[0])
        self.addobjs(o)

        settings["color"] = o[0].get_facecolor()
        
        
class Text(Item):

    def __init__(self,text,p0,p1,**kwargs):
        self._text = text
        self._pdim = [[p0],[p1]]

        super(Text,self).__init__(**kwargs)
    
    @staticmethod
    def defaultsettings():
        return {"color":None,"linewidth":2,"horizontalalignment":"left","verticalalignment":"bottom","xytext":None,"fontsize":15,"fontweight":1}
    
    def datarange(self,dim):
        ran = self.scene.datatransform(self._pdim[dim],dim)
        return min(ran),max(ran)
    
    def pdim(self,dim):
        idim = self.dimindex(dim)
        return self.scene.datatransform(self._pdim[idim],idim)
    
    def update(self):
        scene = self.scene
        settings = scene.getitemsettings(self)
        
        x = scene.xmagnitude(self.pdim(1))[0]
        y = scene.ymagnitude(self.pdim(0))[0]

        if settings["xytext"]!=(0,0):
            arrowprops = dict(facecolor=settings["color"], shrink=0.05)
        else:
            arrowprops = None

        if settings["xytext"] is None:
            xytext = (x,y)
        else:
            xytext = settings["xytext"]
            
        o = [scene.ax.annotate(self._text, xy=(x,y),  xycoords='data',\
                xytext=xytext, textcoords='data',\
                horizontalalignment=settings["horizontalalignment"], verticalalignment=settings["verticalalignment"],\
                color=settings["color"], arrowprops=arrowprops,fontsize=settings["fontsize"],fontweight=settings["fontweight"])]

        self.addobjs(o)

        settings["color"] = o[0].get_color()
        
