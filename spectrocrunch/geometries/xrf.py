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

from ..common.classfactory import with_metaclass
from . import base

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

    def calibrate_fit(self,intensities,\
                detectorposition=None,fit=True,\
                plot=True,xlabel="Motor position (DU)",ylabel="Intensity"):
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
            ax = plt.gcf().gca()
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            plt.legend(loc='best')

        return out
        
        
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

