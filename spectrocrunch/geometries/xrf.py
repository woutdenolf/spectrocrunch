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

from . import base
from ..common.classfactory import with_metaclass
from ..resources import resource_filename

import numpy as np
import silx.math.fit as silxfit
import matplotlib.pyplot as plt


class LinearMotor(object):

    def __init__(self,zerodistance=None,unittocm=None):
        """
        Args:
            zerodistance(num): distance for zero motor position
            unittocm(num): conversion factor motor units -> cm
        """
        self.zerodistance = zerodistance
        self.unittocm = float(unittocm)

    def __call__(self,detectorposition=None,distance=None):
        """Convert detector position to distance and vice versa

        Args:
            detectorposition(Optional(num)): detector position (units)
            distance(Optional(num)): sample-detector distance (cm)
            
        Returns:
            num or None: distance or detector position
        """
        if distance is None:
            return (detectorposition+self.zerodistance)*self.unittocm
        else:
            detectorposition = distance/self.unittocm - self.zerodistance
            return {"detectorposition":detectorposition}

    def calibrate_manual(self,distance,detectorposition=None):
        """Calibrate geometry based on one (distance,position) pair
        """
        self.zerodistance = distance/self.unittocm - detectorposition

    def calibrate_auto(self,rcfile,**kwargs):
        detectorposition,intensities = np.load(self.distance_rcfile)
        return self.calibrate_fit(intensities,detectorposition=detectorposition,**kwargs)

    def calibrate_fit_testcorrelation2(self,intensities,rate=None,zerodistance=None,detectorposition=None):

        # Fit function
        def func(x,rate,zerodistance,activearea):
            distance = (x+zerodistance)*self.unittocm
            sa = self.geometry.solidangle_calc(activearea=activearea,distance=distance)
            return rate*sa

        constraints = [[silxfit.CFIXED,0,0],[silxfit.CFREE,0,0],[silxfit.CFIXED,0,0]]
        
        n = 100
        img = np.zeros((n,n))
        vzerodistance = np.linspace(zerodistance-1,zerodistance+1,n)
        m = 0.9
        vrate = np.linspace(rate*m,rate/m,n)

        for i,zerodistancei in enumerate(vzerodistance):
            for j,ratej in enumerate(vrate):
                obs = func(detectorposition,ratej,zerodistancei,self.geometry.activearea)
                img[i,j] = np.sum((obs-intensities)**2/intensities)

        cax = plt.imshow(img, origin='lower', cmap=plt.cm.jet, interpolation='none', extent=[vrate[0],vrate[-1],vzerodistance[0],vzerodistance[-1]])
        ax = plt.gcf().gca()
        ax.set_ylabel("$x_0$ (cm)")
        ax.set_xlabel("$c_x$ (sr$^{-1}$)")
        ax.set_aspect(abs(vrate[-1]-vrate[0])/abs(vzerodistance[-1]-vzerodistance[0]))
        ax.axvline(x=rate)
        ax.axhline(y=zerodistance)
        cbar = plt.colorbar(cax,label="$\chi^2$")
        ax = plt.gcf().gca()
        
    def calibrate_fit_testcorrelation(self,intensities,rate=None,zerodistance=None,detectorposition=None):

        # Fit function
        def func(x,rate,zerodistance,activearea):
            distance = (x+zerodistance)*self.unittocm
            sa = self.geometry.solidangle_calc(activearea=activearea,distance=distance)
            return rate*sa

        constraints = [[silxfit.CFIXED,0,0],[silxfit.CFREE,0,0],[silxfit.CFIXED,0,0]]
        
        n = 50
        y = np.zeros(n)
        y2 = np.zeros(n)
        x = np.linspace(self.geometry.activearea-0.01,self.geometry.activearea+0.01,n)
        for i,activearea in enumerate(x):
            p0 = [rate,self.zerodistance+np.random.uniform(-1,1),activearea]
            p, cov_matrix, info = silxfit.leastsq(func, detectorposition, intensities, p0,\
                                              constraints=constraints, full_output=True)
            y[i] = info["reduced_chisq"]
            y2[i] = p[1]
            
        p = plt.plot(x*100,y)
        ax = plt.gcf().gca()
        color = p[-1].get_color()
        ax.axvline(x=self.geometry.activearea*100,linestyle='dashed',color=color)
        ax.set_xlabel("Active area ($mm^2$)")
        ax.set_ylabel("Reduced-$\chi^2$",color=color)
        ax.tick_params(axis='y', labelcolor=color)
        
        color = next(ax._get_lines.prop_cycler)['color']
        ax2 = ax.twinx()
        ax2.plot(x*100,y2*10,color=color)
        ax2.axhline(y=zerodistance*10,linestyle='dashed',color=color)
        ax2.set_ylabel("$x_0$ (mm)",color=color)
        ax2.tick_params(axis='y', labelcolor=color)
         
    def calibrate_fit(self,intensities,\
                detectorposition=None,fit=True,\
                plot=False,xlabel="Motor position (DU)",ylabel="Intensity"):
        """Calibrate geometry based in intensity vs. linear motor position
        """
        
        # Initial parameter values
        activearea = self.geometry.activearea
        zerodistance = self.zerodistance
        
        distance = self(detectorposition=detectorposition)
        rate = np.mean(intensities/self.geometry.solidangle_calc(activearea=activearea,distance=distance))
        
        p0 = [rate,zerodistance,activearea]
        constraints = [[silxfit.CFREE,0,0],[silxfit.CFREE,0,0],[silxfit.CFIXED,0,0]]

        # Fit function
        def func(x,rate,zerodistance,activearea):
            distance = (x+zerodistance)*self.unittocm
            sa = self.geometry.solidangle_calc(activearea=activearea,distance=distance)
            return rate*sa

        if fit:
            # Fit
            p, cov_matrix, info = silxfit.leastsq(func, detectorposition, intensities, p0,\
                                              constraints=constraints, full_output=True)

            # Parse result                       
            errors = np.sqrt(np.diag(cov_matrix))
            S = np.diag(1/errors)
            cor_matrix = S.dot(cov_matrix).dot(S)

            out =  "rate = {} +/- {}\n".format(p[0],errors[0])
            out += "zerodistance = {} +/- {} DU\n".format(p[1],errors[1])
            out += "activearea = {} +/- {} cm^2\n".format(p[2],errors[2])
            out += "R(rate,zerodistance) = {}\n".format(cor_matrix[0,1])
            out += "R(rate,activearea) = {}\n".format(cor_matrix[0,2])
            out += "R(zerodistance,activearea) = {}\n".format(cor_matrix[1,2])
        
            rate,zerodistance,activearea = p
        
            # Apply result
            self.zerodistance = zerodistance
            #self.geometry.activearea = activearea
            
            label = 'fit (active area fixed)'
        else:
            out = ""
            label = 'current geometry'

        # Plot
        if plot:
            plt.plot(detectorposition,intensities,'x',label='data')
            plt.plot(detectorposition,func(detectorposition,rate,zerodistance,activearea),label=label)
            off = 0.6
            plt.figtext(0.55,off,"$x_0$ = {} mm".format(self.zerodistance*10))
            plt.figtext(0.55,off-0.05,"$A$ = {} $mm^2$".format(self.geometry.activearea*100))
            plt.figtext(0.55,off-0.1,"$c_x$ = {}".format(rate))
            
            ax = plt.gcf().gca()
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            plt.legend(loc='best')

        return {"rate":rate,"zerodistance":zerodistance,"activearea":activearea},out
        
        
class Geometry(with_metaclass(base.Point)):

    def __init__(self,distancefunc=None,distanceargs=None,**kwargs):
        """
        Args:
            distancefunc(Optional(callable)): distance = distancefunc(**distanceargs)
            distanceargs(Optional(dict)): 
        """

        self.set_distancefunc(distancefunc)
        self.distanceargs = distanceargs
        
        super(Geometry,self).__init__(**kwargs)

    def set_distancefunc(self,distancefunc):
        if distancefunc is not None:
            distancefunc.geometry = self
        self.distancefunc = distancefunc

    @property
    def distance(self):
        """Sample-detector distance in cm
        """
        if self.distancefunc is None:
            return self._distance
        else:
            return self.distancefunc(**self.distanceargs)
        
    @distance.setter
    def distance(self,distance):
        if self.distancefunc is None:
            self._distance = distance
        else:
            self.distanceargs.update(self.distancefunc(distance=distance))
    
    def calibrate_distance_manual(self,distance,**distanceargs):
        if self.distancefunc is not None:
            self.distancefunc.calibrate_manual(distance,**distanceargs)

    def calibrate_distance_fit(self,intensities,**kwargs):
        if self.distancefunc is not None:
            return self.distancefunc.calibrate_fit(intensities,**kwargs)
    
    @property
    def distance_rcfile(self):
        return resource_filename('geometry/distancecalib_{}_{}.npy'.format(self.__class__.__name__,self.detector.__class__.__name__))
    
    def calibrate_distance_auto(self,**kwargs):
        if self.distancefunc is not None:
            return self.distancefunc.calibrate_auto(self.distance_rcfile,**kwargs)
        
    
class sxm120(Geometry):

    def __init__(self,**kwargs):
        # blc10516 (April 2017)
        # detector position in mm
        
        zerodistance = kwargs.pop("zerodistance",60.38)
        distancefunc = LinearMotor(zerodistance=zerodistance,unittocm=0.1)

        super(sxm120,self).__init__(anglein=62,angleout=49,azimuth=0,\
                        distancefunc=distancefunc,**kwargs)

class sxm90(Geometry):

    def __init__(self,**kwargs):
        # blc10516 (April 2017)
        # detector position in mm
        
        zerodistance = kwargs.pop("zerodistance",85.49)
        distancefunc = LinearMotor(zerodistance=zerodistance,unittocm=0.1)

        super(sxm90,self).__init__(anglein=62,angleout=28,azimuth=0,\
                        distancefunc=distancefunc,**kwargs)

class microdiff(Geometry):

    def __init__(self,**kwargs):
        # blc10516 (April 2017)
        # detector position in mm
        
        zerodistance = kwargs.pop("zerodistance",-57.30)
        distancefunc = LinearMotor(zerodistance=zerodistance,unittocm=-0.1)

        super(microdiff,self).__init__(anglein=90,azimuth=0,\
                        distancefunc=distancefunc,**kwargs)
                        

factory = Geometry.factory

