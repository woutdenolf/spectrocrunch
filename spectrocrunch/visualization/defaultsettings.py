# -*- coding: utf-8 -*-
#
#   Copyright (C) 2018 European Synchrotron Radiation Facility, Grenoble, France
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
import matplotlib
import logging

from ..common import instance

logger = logging.getLogger(__name__)

def update(parameters):
    matplotlib.rcParams.update(parameters)

def figsize(publish="screen",aspect=0.75,nsidebyside=1,space=0.,widescreen=True):
    """
    Args:
        publish(Optional(str)): word, powerpoint, screen
        aspect(Optional(num)): V/H
        nsidebyside(Optional(num)): number of images that should be put side-by-side
        space(Optional(num)): fraction of the total width that is white space
        widescreen(Optional(bool)): only for powerpoint
    """
    if publish=="word":
        width = 6.24
    elif publish=="powerpoint":
        if widescreen:
            width = 33.867
        else:
            width = 13.33
    else:
        if instance.isnumber(publish):
            width = publish
        else:
            width = 6.4
    width *= (1-space)/nsidebyside
    return width,width*aspect
        
def adapttofigsize(size,fontsize=None,**kwargs):
    """
    Args:
        size(2-tuple): inch
        fontsize(Optional(num)): inch
    """
    if fontsize: # in inch
        # 1 inch = 72 points
        fontsize = fontsize*72
    else:
        fontsize = min(size)*10/4.8
    linewidth = fontsize*0.08
    
    wtick = linewidth
    cwtick = 0.75
    
    stick = linewidth*4.375
    cstick = 1/1.75
    
    ptick = linewidth*4.375
    cptick = 3.4/3.5
    
    parameters={"font.size":fontsize,
              "axes.linewidth":linewidth,
              "lines.linewidth":linewidth*1.875,
              "patch.linewidth":linewidth*1.25,
              "axes.titlepad":fontsize*0.6,
              "axes.labelpad":fontsize*0.4,
              "xtick.major.width":wtick,
              "xtick.minor.width":wtick*cwtick,
              "ytick.major.width":wtick,
              "ytick.minor.width":wtick*cwtick,
              "xtick.major.size":stick,
              "xtick.minor.size":stick*cstick,
              "ytick.major.size":stick,
              "ytick.minor.size":stick*cstick,
              "xtick.major.pad":ptick,
              "xtick.minor.pad":ptick*cptick,
              "ytick.major.pad":ptick,
              "ytick.minor.pad":ptick*cptick,
              "xtick.direction":kwargs.get("tick.direction","out"),
              "ytick.direction":kwargs.get("tick.direction","out"),
                }
    
    update(parameters)
    
    parameters = {k:v for k,v in kwargs.items() if k in matplotlib.rcParams}
    update(parameters)

def screensize():
    figs1 = plt.get_fignums()
    mgr = plt.get_current_fig_manager()
    h = mgr.window.winfo_screenheight()
    w = mgr.window.winfo_screenwidth()
    figs2 = plt.get_fignums()
    if len(figs2)>len(figs1):
        plt.close(figs2[-1])
    return w,h

def dpi(publish="photo&text",best=True):
    if publish=="color":
        ret = 300
    elif publish=="b&w":
        ret = 600 if best else 300
    elif publish=="photo&text":
        ret = 900 if best else 600
    elif publish=="lineart":
        ret = 1200 if best else 900
    else: # width x height
        if not publish or not instance.isarray(publish):
            publish = screensize()
        width,height = publish
        
        # Powerpoint: normal and widescreen
        #50	 500 × 375   667 × 375   50 dpi
        #96  960 × 720   1280 × 720  96 dpi
        #100 1000 × 750  1333 × 750  100 dpi
        #150 1500 × 1125 2000 × 1125 150 dpi
        #200 2000 × 1500 2667 × 1500 200 dpi
        #250 2500 × 1875 3333 × 1875 250 dpi
        #300 3000 × 2250 4000 × 2250 300 dpi
        if best:
            warr = [667,1280,133,200,2667,333,4000]
        else:
            warr = [500,960,1000,1500,2000,2500,3000]
        harr = [375,720,750,1125,1500,1875,2250]
        dpiarr = [50,96,100,150,200,250,300]
        
        for w,h,r in zip(warr,harr,dpiarr):
            if w>width and h>height:
                ret = r
                break
        else:
            ret = 300

    logger.info("DPI {} chosen to publish \"{}\"".format(ret,publish))

    return ret

