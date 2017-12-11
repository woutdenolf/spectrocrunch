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
import numpy as np
from ..common import instance
from ..math import fit1d

def plotline(energy,intensity,FWHM=None,out=False,name=None,**kwargs):
    if FWHM is None:
        plt.plot([energy,energy],[0,intensity],**kwargs)
        h = intensity
    else:
        sx = FWHM/(2*np.sqrt(2*np.log(2)))
        k = 4
        x = np.linspace(energy-k*sx,energy+k*sx,50)
        y = fit1d.gaussian(x,energy,sx,intensity)
        h = max(y)
        plt.plot(x,y,**kwargs)
    
    if name is not None:
        kwargs2 = {}
        if "color" in kwargs:
            kwargs2["color"] = kwargs["color"]
        
        plt.annotate(name, xy=(energy, h),**kwargs2)

def xrf(spectrum,FWHM=None,out=False,mark=True):
    ax = plt.gca()
    
    for element,shells in spectrum["fluorescence"].items():
        for shell,lines in shells.items():
            label = "{} {}".format(element,shell)
            color = next(ax._get_lines.prop_cycler)['color']
            if mark:
                mark = lines.keys()[np.argmax(lines.values())]
            for line,intensity in lines.items():
                energy = line.energy(element)
                if out:
                    print element,line,energy,intensity
                if mark==line:
                    name = "{} {}".format(element,line)
                else:
                    name = None
                plotline(energy,intensity,FWHM=FWHM,name=name,color=color,label=label)
                label = None
    
    scatangle = 111
    
    for line,intensity in spectrum["scattering"].items():
        label = str(line)
        color = next(ax._get_lines.prop_cycler)['color']
        energy = line.energy(scatangle)
        if out:
            print line,energy,intensity
        if instance.isarray(intensity):
            if mark:
                mark = lines.keys()[np.argmax(lines.values())]
            for energy,intensity in zip(energy,intensity):
                print mark
                if mark==line:
                    name = str(line)
                else:
                    name = None
                plotline(energy,intensity,FWHM=FWHM,color=color,name=name,label=label)
                label = None
        else:
            if mark:
                name = str(line)
            else:
                name = None
            plotline(energy,intensity,FWHM=FWHM,name=name,color=color,label=label)
    
    #if FWHM is not None:
    #    ax = plt.gca()
    #    ax.set_yscale('log', basey=10)
    #    plt.ylim([1,None])
        
    plt.legend(loc='best')
    plt.xlabel("Energy (keV)")
    if "ylabel" in spectrum:
        plt.ylabel(spectrum["ylabel"])
    if "xlim" in spectrum:
        plt.xlim(spectrum["xlim"])
    if "title" in spectrum:
        plt.title(spectrum["title"])

