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
from ..common import units
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
    
    
class Geometry2D(object):
    
    def dataaxis(self,displayaxis):
        """
        Args:
            displayaxis(num): display axis (0 or 1)
            
        Returns:
            dataaxis(num): data axis (0 or 1)
        """
        if instance.isstring(displayaxis):
            displayaxis = int(displayaxis=='x')
            
        if self.transposed:
            return 1-displayaxis
        else:
            return displayaxis

    @property
    def dataaxisx(self):
        """Data axis corresponding to the X-axis
        """
        return self.dataaxis('x')

    @property
    def dataaxisy(self):
        """Data axis corresponding to the Y-axis
        """
        return self.dataaxis('y')


class Scene(Geometry2D):
    """Each scene can have a number of registered items with associalted plot settings.
    """

    def __init__(self,ax=None,unit0=ureg.dimensionless,unit1=ureg.dimensionless,title=""):
        self.ax = ax
        self._items = OrderedDict()
        self.transposed = False
        self.cmap = plt.get_cmap('jet')
        self.title = title
        self._aspectratio_display = None
        
        # Order: data dimensions
        self.axlabels = ["Dim0","Dim1"]
        self._increasing = [True,True]
        self.units = [ureg.Unit(unit0),ureg.Unit(unit1)]
        self.dataoffset = [ureg.Quantity(0,unit0),ureg.Quantity(0,unit1)]
        self.datascale = [1,1]
    
    def q0(self,v):
        return self.dimquantity(v,0)
    
    def q1(self,v):
        return self.dimquantity(v,1)
          
    def dimquantity(self,v,dataaxis):
        if isinstance(v,ureg.Quantity):
            return v
        else:
            return ureg.Quantity(v,self.units[dataaxis])

    def xmagnitude(self,v):
        return self._magnitude(v,self.dataaxisx)
    
    def ymagnitude(self,v):
        return self._magnitude(v,self.dataaxisy)
        
    def _magnitude(self,v,dataaxis):
        unit = self.units[dataaxis]
        if instance.isarray(v):
            return type(v)(a.to(unit).magnitude for a in v)
        else:
            return v.to(unit).magnitude

    @property
    def xunit(self):
        return self.units[self.dataaxisx]
    
    @xunit.setter
    def xunit(self,value):
        self.units[self.dataaxisx] = ureg.Unit(value)

    @property
    def yunit(self):
        return self.units[self.dataaxisy]
    
    @yunit.setter
    def yunit(self,value):
        self.units[self.dataaxisy] = ureg.Unit(value)
        
    @property
    def xlabel(self):
        return self._label(self.dataaxisx)
    
    @property
    def ylabel(self):
        return self._label(self.dataaxisy)
        
    def _label(self,dataaxis):
        unit = self.units[dataaxis]
        if unit == ureg.dimensionless:
            return self.axlabels[dataaxis]
        else:
            unit = "{:~}".format(unit)
            if unit=="um":
                unit="$\mu$m"
            return "{} ({})".format(self.axlabels[dataaxis],unit)

    def setaxes(self,ax):
        self.ax = ax
    
    def showaxes(self):
        self.ax.set_axis_on()
        self.ax.get_dataaxisx().set_visible(True)
        self.ax.get_dataaxisy().set_visible(True)
    
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
    
    def __iter__(self):
        for item in self._items:
            item.selectscene(self)
        return iter(self._items)
    
    @property
    def items(self):
        return list(self)
    
    def datarange(self,dataaxis):
        ran = [r for item in self for r in item.datarange(dataaxis)]
        return min(ran),max(ran)

    def setdatarange(self,dataaxis,ran):
        self.dataoffset[dataaxis] = self.dimquantity(0,dataaxis)
        self.datascale[dataaxis] = 1
        ranorg = sorted(self.datarange(dataaxis))
        scale = (ran[1]-ran[0])/(ranorg[1]-ranorg[0])
        off = scale*ranorg[0]-ran[0]
        self.dataoffset[dataaxis] = off
        self.datascale[dataaxis] = scale
    
    def resetdatarange(self):
        self.dataoffset = [self.q0(0),self.q1(1)]
        self.datascale = [1,1]

    def increasingx(self,increasing):
        self._increasing[self.dataaxisx] = increasing

    def increasingy(self,increasing):
        self._increasing[self.dataaxisy] = increasing
        
    def flipx(self,increasing=None):
        self._flipaxis(self.dataaxisx,increasing=increasing) 

    def flipy(self,increasing=None):
        self._flipaxis(self.dataaxisy,increasing=increasing) 
        
    def _flipaxis(self,dataaxis,increasing=None):
        if increasing is None:
            self._increasing[dataaxis] = not self._increasing[dataaxis]
        else:
            if increasing==self._increasing[dataaxis]:
                self.setdatarange(dataaxis,self.datarange(dataaxis)[::-1])
            else:
                self._increasing[dataaxis] = not self._increasing[dataaxis]

    def origin(self,dataaxis,off):
        self.dataoffset[dataaxis] = off-self.dataoffset[dataaxis]
    
    def originreset(self):
        self.origin(0,self.q0(0))
        self.origin(1,self.q1(0))
    
    def originzero(self):
        self.origin(0,self.datarange(0)[0])
        self.origin(1,self.datarange(1)[0])
        
    def datatransform(self,arr,dataaxis):
        arr,func = units.asarrayf(arr)
        return func(arr*self.datascale[dataaxis] - self.dataoffset[dataaxis])
    
    def displaylim(self,lim,dataaxis):
        if self._increasing[dataaxis]:
            return min(lim),max(lim)
        else:
            return max(lim),min(lim)
    
    @property
    def displaylimx(self):
        datalimx = [x for item in self for x in item.datalimx]
        return self.displaylim(datalimx,self.dataaxisx)

    @property
    def displaylimy(self):
        ylim = [y for item in self for y in item.datalimy]
        return self.displaylim(ylim,self.dataaxisy)

    @property
    def aspectcorrect(self):
        if self._aspectratio_display is None:
            return None # same as 1: aspect ratio of display and data are the same
        return self._aspectratio_display/self.aspectratio_data
    
    @property
    def aspectratio(self):
        if self._aspectratio_display is None:
            return self.aspectratio_data
        else:
            return self._aspectratio_display
            
    @aspectratio.setter
    def aspectratio(self,value):
        self._aspectratio_display = value
        
    @property
    def aspectratio_data(self):
        displaylimy = self.displaylimy
        displaylimx = self.displaylimx
        dx = abs(displaylimx[1]-displaylimx[0])
        dy = abs(displaylimy[1]-displaylimy[0])
        return (dy/dx).to("dimensionless").magnitude
    
    def crop(self):
        lim = self.xmagnitude(self.displaylimx)
        self.ax.set_xlim(lim[0],lim[1])
        #self.ax.set_xbound(lim[0],lim[1])
        lim = self.ymagnitude(self.displaylimy)
        self.ax.set_ylim(lim[0],lim[1])
        #self.ax.set_ybound(lim[0],lim[1])
        
    def update(self,**kwargs):
        for item in self:
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
        
        tmp = [single(item.vminmax) for item in self if hasattr(item,'vminmax')]
        
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
    
    @property
    def wsize_pixels(self):
        pts = self.ax.get_window_extent().get_points()
        x1,x2 = pts[:,1]
        y1,y2 = pts[:,1]
        return x2-x1+1,y2-y1+1

    @property
    def wsize_data(self):
        x1,x2 = self.xmagnitude(self.displaylimx)
        y1,y2 = self.ymagnitude(self.displaylimy)
        return x2-x1,y2-y1
    
    def fontsize_data(self,pt):
        xp,yp = self.wsize_pixels
        xd,yd = self.wsize_data
        return pt*xd/xp,pt*yd/yp
           
           
class Item(Hashable,Geometry2D):
    """Each item can be registered with multiple scenes. Data is owned by the Item instance,
       plot settings are owned by the scenes.
    """
    
    def __init__(self,axis0name="Dim0",axis1name="Dim1",scene=None,name=None,**kwargs):
        self._scenes = []
        self._plotobjs = [] # list of lists, each list belonging to a scene
        self._sceneindex = -1
        self.axis0name = axis0name
        self.axis1name = axis1name
        self.name = name
        
        if scene is not None:
            self.register(scene)
        
        try:
            self.set_settings(kwargs)
        except RuntimeError:
            logger.warning("Item settings are not applied (provide a scene)")
    
    def register(self,scene):
        scene.register(self)
    
    def useaxesnames(self):
        self.scene.axlabels = [self.axis0name,self.axis1name]
    
    def defaultsettings(self):
        return {}
    
    def set_settings(self,settings):
        settings2 = self.scene.getitemsettings(self)
        for k,v in settings.items():
            settings2[k] = v

    def get_settings(self, keys):
        settings = self.scene.getitemsettings(self)
        return {k:settings[k] for k in keys}
        
    def set_setting(self,attr,value):
        settings = self.scene.getitemsettings(self)
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
            

    def datarange(self,dataaxis):
        raise NotImplementedError("Item is an abstract class")
    
    @property
    def datalimx(self):
        return self.datarange(self.dataaxisx)
    
    @property
    def datalimy(self):
        return self.datarange(self.dataaxisy)
    
    def _stringrepr(self):
        return "{} {}".format(type(self).__name__,id(self))
    
    
class Image(Item):

    def __init__(self,data,lim0=None,lim1=None,labels=None,**kwargs):
        self.data = data
        self.labels = labels
        
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
    def data(self):
        """Data for the selected channels
        """
        channels = self.channels
        shape = (self._data.shape[0],self._data.shape[1],len(channels))
        data = np.zeros(shape,self.dtype)
        for ch,i in enumerate(channels):
            if i is not None:
                data[...,ch] = self._get_data(i)
        if len(channels)==1:
            data = data[...,0]
        return data

    def _get_data(self,index):
        if self._data.ndim==3:
            return self._data[...,index]
        else:
            if index!=0:
                raise IndexError("Raw data has only one channel, you tried to select channel {}".format(index))
            return self._data
            
    @data.setter
    def data(self,value):   
        self._data = value

    @property
    def dtype(self):
        return self._data.dtype

    @property
    def labels(self):
        """Labels for the selected channels
        """
        labels = [None]*self.nchannels
        for ch,i in enumerate(self.channels):
            if i is not None:
                try:
                    labels[ch] = self._labels[i]
                except:
                    labels[ch] = "<empty>"
        return labels

    @labels.setter
    def labels(self,value):
        self._labels = value

    @property
    def channels(self):
        """Return 1 or 3 channels
        """
        settings = self.scene.getitemsettings(self)
        channels = settings["channels"]
        if channels is None:
            nchannels = min(self.datadepth,3)
            channels = range(nchannels)
        else:
            if not instance.isarray(channels):
                channels = [channels]
            nchannels = min(len(channels),3)
        if nchannels!=1:
            nchannels=3
        
        if len(channels)<nchannels:
            channels.extend([None]*(nchannels-len(channels)))
        else:
            channels = channels[0:nchannels]
        return channels
    
    @property
    def datadepth(self):
        if self._data.ndim==3:
            return self._data.shape[-1]
        else:
            return 1
            
    @property
    def nchannels(self):
        return len(self.channels)
    
    @property
    def datadisplay(self):
        if self.transposed:
            return np.swapaxes(self.data,0,1)
        else:
            return self.data

    @property
    def vminmax(self):
        img = self.data
        if self.data.ndim==3:
            nchannels = self.data.shape[-1]
            vmin = np.asarray([np.nanmin(img[...,i]) for i in range(nchannels)])
            vmax = np.asarray([np.nanmax(img[...,i]) for i in range(nchannels)])
        else:
            vmin = np.nanmin(img)
            vmax = np.nanmax(img)
        return vmin,vmax
    
    def scale(self,vmin=None,vmax=None):
        self.set_setting("vmin",vmin)
        self.set_setting("vmax",vmax)
    
    @staticmethod
    def _pminmaxparse(p,v):
        # Handle mismatch in number of channels
        p = instance.asscalar(p)
        if instance.isarray(p):
            p = instance.asarray(p)
            if instance.isarray(v):
                if p.size<v.size:
                    p = np.append(p,[0]*(v.size-p.size))
                else:
                    p = p[0:v.size]
            else:
                p = p[0]
        return p,v
        
    def selfscale(self,pmin=0,pmax=1):
        vmin,vmax = self.vminmax
        d = vmax-vmin
        pmin,vmin = self._pminmaxparse(pmin,vmin)
        pmax,vmax = self._pminmaxparse(pmax,vmax)
        vmax = vmin + d*pmax
        vmin = vmin + d*pmin
        self.scale(vmin=vmin,vmax=vmax)
    
    @staticmethod
    def defaultsettings():
        return {"cmap":None,\
                  "aspectcorrect":None,\
                  "border":False,\
                  "color":None,\
                  "image":True,\
                  "linewidth":2,\
                  "alpha":1,\
                  "vmin":None,\
                  "vmax":None,\
                  "cnorm":None,\
                  "legend":True,\
                  "fontsize":12,\
                  "fontweight":500,\
                  "legendposition":"RT",\
                  "channels":None}
        
    def datarange(self,dataaxis,border=False):
        lim = self.scene.datatransform(self.lim[dataaxis],dataaxis)
        if border:
            d = (lim[1]-lim[0])/(2.*(self._data.shape[dataaxis]-1))
            return lim[0]-d,lim[1]+d
        else:
            return lim[0],lim[1]
    
    @property
    def datalimx(self):
        return self.datarange(self.dataaxisx,border=True)
    
    @property
    def datalimy(self):
        return self.datarange(self.dataaxisy,border=True)

    @property
    def displaylimx(self):
        scene = self.scene
        return scene.displaylim(self.datalimx,scene.dataaxisx)

    @property
    def displaylimy(self):
        scene = self.scene
        return scene.displaylim(self.datalimy,scene.dataaxisy)
    
    @property
    def labelxy(self):
        scene = self.scene
        settings = scene.getitemsettings(self)
        
        if "S" in settings["legendposition"]:
            o = scene
        else:
            o = self
            
        if "L" in settings["legendposition"]:
            x = 0
            horizontalalignment = "right"
            ox = -0.5
            dx = 0
        else:
            x = 1
            horizontalalignment = "left"
            ox = 0.5
            dx = 0
            
        if "B" in settings["legendposition"]:
            y = 0
            verticalalignment = "bottom"
            oy = 0
            dy = 1.5
        else:
            y = 1
            verticalalignment = "top"
            oy = 0
            dy = -1.5

        x = scene.xmagnitude(o.displaylimx[x])
        y = scene.ymagnitude(o.displaylimy[y])
        xchar,ychar = scene.fontsize_data(settings["fontsize"])

        return x+xchar*ox,y+ychar*oy,xchar*dx,ychar*dy,horizontalalignment,verticalalignment

    def bordercoord(self):
        x = sorted(self.datalimy)
        y = sorted(self.datalimy)
        x = [x[0],x[1],x[1],x[0],x[0]]
        y = [y[0],y[0],y[1],y[1],y[0]]
        return x,y

    def update(self):
        scene = self.scene
        
        settings = scene.getitemsettings(self)
        if settings["cmap"] is None:
            cmap = scene.cmap
        else:
            cmap = settings["cmap"]
        if settings["aspectcorrect"] is None:
            aspectcorrect = scene.aspectcorrect
        else:
            aspectcorrect = settings["aspectcorrect"]
        alpha = settings["alpha"]
        
        # Create intensity-scaled image
        vmin,vmax = scene.vminmax
        if settings["vmin"] is not None:
            vmin = settings["vmin"]   
        if settings["vmax"] is not None:
            vmax = settings["vmax"]
            
        image = self.datadisplay
        
        if image.ndim==3:
            nchannels = image.shape[-1]
            if not instance.isarray(vmin):
                vmin = [vmin]*nchannels
            if not instance.isarray(vmax):
                vmax = [vmax]*nchannels
            for i in range(nchannels):
                if settings["cnorm"] is None:
                    norm = pltcolors.Normalize(vmin=vmin[i],vmax=vmax[i],clip=True)
                else:
                    norm = settings["cnorm"](vmin=vmin[i],vmax=vmax[i],clip=True)
                    
                image[...,i] = norm(image[...,i]).data
            norm = None
        else:
            nchannels = 1
            if settings["cnorm"] is None:
                norm = pltcolors.Normalize(vmin=vmin,vmax=vmax)
            else:
                norm = settings["cnorm"](vmin=vmin,vmax=vmax)
        
        # Coordinates and borders
        extent = scene.xmagnitude(self.datalimx)+scene.ymagnitude(self.datalimy)
        
        x,y = self.bordercoord()
        x = scene.xmagnitude(x)
        y = scene.ymagnitude(y)
        
        # Legend (right-top corner)
        xlabel,ylabel,dx,dy,horizontalalignment,verticalalignment = self.labelxy

        # Create/update objects
        o = self.plotobjs
        update = bool(o)
        
        # Image + border
        if update:
            o[0].set_data(image)
            plt.setp(o[0], interpolation = 'nearest',\
                        extent = extent,\
                        norm = norm,\
                        alpha = alpha,\
                        cmap = cmap)
            o[0].set_visible(settings["image"])
            
            o[1].set_data(x,y)
            kwargs = {k:settings[k] for k in ["linewidth","color","alpha"]}
            plt.setp(o[1], **kwargs)
            o[1].set_visible(settings["border"])
            settings["color"] = o[1].get_color()
            
        else:
            # matplotlib.image.AxesImage
            self.cax = oi = scene.ax.imshow(image,\
                        interpolation = 'nearest',\
                        extent = extent,\
                        cmap = cmap,\
                        origin = 'lower',
                        norm = norm,\
                        alpha = alpha,\
                        aspect = aspectcorrect)
            oi.set_visible(settings["image"])
            o = [oi]
            
            kwargs = {k:settings[k] for k in ["linewidth","color","alpha"]}
            oi = scene.ax.plot(x,y,**kwargs)[0]
            oi.set_visible(settings["border"])
            settings["color"] = oi.get_color()
            o.append(oi)
        
        # Legend
        kwargs = {k:settings[k] for k in ["color","alpha"]}
        kwargs["horizontalalignment"] = horizontalalignment
        kwargs["verticalalignment"] = verticalalignment

        if nchannels==1:
            colors = [kwargs["color"]]
        else:
            colors = ['r','g','b']
            if nchannels>3:
                colors.extend([None]*(nchannels-3))
            colors = colors[0:nchannels]

        i = -1
        for label,color in zip(self.labels,colors):
            if label is None or color is None:
                continue
            else:
                i += 1
            kwargs["color"] = color
            if update:
                oindex = i+2
                o[oindex].set_text(label)
                o[oindex].set_position((xlabel+i*dx,ylabel+i*dy))
                plt.setp(o[oindex], **kwargs)
                o[oindex].set_visible(settings["legend"])
            else:
                oi = scene.ax.text(xlabel+i*dx,ylabel+i*dy,label,**kwargs)
                oi.set_visible(settings["legend"])
                o.append(oi)
                
        if not update:
            self.addobjs(o)

            
class Scatter(Item):

    def __init__(self,coordinate0,coordinate1,labels=None,**kwargs):
        self._coordinate = [coordinate0,coordinate1]
        if labels is None:
            self.labels = []
        else:
            self.labels = labels
        super(Scatter,self).__init__(**kwargs)
    
    def datarange(self,dataaxis):
        ran = self.scene.datatransform(self._coordinate[dataaxis],dataaxis)
        return min(ran),max(ran)
    
    def coordinate(self,dataaxis):
        return self.scene.datatransform(self._coordinate[dataaxis],dataaxis)
        
    @staticmethod
    def defaultsettings():
        return {"color":None,"marker":"+","linestyle":"","linewidth":2,\
                "fill":False,"alpha":None,"closed":False,"labels":True,\
                "horizontalalignment":"left","verticalalignment":"bottom",\
                "fontsize":12,"labeloffset":0.1,"fontweight":500}

    def update(self):
        scene = self.scene
        settings = scene.getitemsettings(self)
        
        x = instance.asarray(scene.xmagnitude(self.coordinate(self.dataaxisx)))
        y = instance.asarray(scene.ymagnitude(self.coordinate(self.dataaxisy)))

        o = self.plotobjs

        if settings["fill"]:
            update = len(o)==1
            if update:
                update = isinstance(o[0],matplotlib.patches.Polygon)
        
            kwargs = {k:settings[k] for k in ["linestyle","linewidth","color","alpha","closed"]}
            if not kwargs["linestyle"]:
                kwargs["linestyle"] = None
            
            # Create/update polygon
            if update:
                o[0].set_xy=zip(x, y)
                plt.setp(o[0],**kwargs)
            else:
                o = [matplotlib.patches.Polygon(zip(x, y),**kwargs)]
                scene.ax.add_patch(o[0])

            settings["color"] = o[0].get_facecolor()
        
        elif settings["linestyle"] or settings["marker"]:
            blabels = settings["labels"] and bool(self.labels)
            update = len(o) == len(self.labels)*blabels+1
            
            if settings["closed"] and settings["linestyle"]:
                xcl = np.append(x,x[0])
                ycl = np.append(y,y[0])
            else:
                xcl = x
                ycl = y

            # Create/update markers and/or lines
            kwargs = {k:settings[k] for k in ["marker","linestyle","alpha","color","linewidth"]}
            if update:
                o[0].set_data(xcl,ycl)
                plt.setp(o[0],**kwargs)
                settings["color"] = o[0].get_color()
            else:
                o = [scene.ax.plot(xcl,ycl,**kwargs)[0]]
                settings["color"] = o[0].get_color()

            # Create labels
            if blabels:
                xo,yo = scene.fontsize_data(settings["fontsize"])
                xo *= settings["labeloffset"]
                yo *= settings["labeloffset"]

                kwargs = {k:settings[k] for k in ["horizontalalignment","verticalalignment","alpha","color","fontsize","fontweight"]}
                if not update:
                    kwargs["xycoords"] = "data"
                    kwargs["textcoords"] = "data"
                
                for i,(xi,yi,label) in enumerate(zip(x,y,self.labels)):
                    if update:
                        o[i+1].set_text(label)
                        o[i+1].set_position((xi,yi))
                        o[i+1].xytext = (xi+xo,yi+yo) # has no effect!!!!
                        plt.setp(o[i+1],**kwargs)
                    else:
                        kwargs["xy"] = (xi,yi)
                        kwargs["xytext"] = (xi+xo,yi+yo)
                        o.append(scene.ax.annotate(label,**kwargs))
        if not update:
            self.addobjs(o)


class Polyline(Scatter):

    @staticmethod
    def defaultsettings():
        settings = Scatter.defaultsettings()
        settings["marker"] = ""
        settings["linestyle"] = "-"
        settings["closed"] = True
        return settings


class Polygon(Scatter):

    @staticmethod
    def defaultsettings():
        settings = Scatter.defaultsettings()
        settings["marker"] = ""
        settings["fill"] = True
        settings["closed"] = True
        return settings
        
        
class Text(Scatter):
        
    @staticmethod
    def defaultsettings():
        settings = Scatter.defaultsettings()
        settings["labels"] = True
        return settings

