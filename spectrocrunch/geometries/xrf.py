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
from ..simulation import noisepropagation
from ..common import units

import numpy as np
import silx.math.fit as silxfit
import matplotlib.pyplot as plt
import os
import json

class LinearMotor(object):

    def __init__(self,zerodistance=None,positionsign=1):
        self.zerodistance = zerodistance
        self.positionsign = positionsign

    @property
    def zerodistance_rv(self):
        return self._zerodistance
    
    @property
    def zerodistance(self):
        return noisepropagation.E(self.zerodistance_rv)
    
    @zerodistance.setter
    def zerodistance(self,value):
        if value is None:
            self._zerodistance = None
        else:
            self._zerodistance = units.Quantity(value,"cm")

    def _distancecalc(self,detectorposition=None,zerodistance=None,distance=None):
        if distance is None:
            return self.positionsign * (detectorposition+zerodistance)
        elif detectorposition is None:
            distance = units.Quantity(distance,"cm")
            return self.positionsign*distance - zerodistance
        elif zerodistance is None:
            distance = units.Quantity(distance,"cm")
            return self.positionsign*distance - self.geometry.detectorposition
        else:
            raise RuntimeError("Either distance, detector position or zero-distance should be unknown")
        
    def __call__(self,detectorposition=None,distance=None):
        """Convert detector position to distance and vice versa

        Args:
            detectorposition(Optional(num)): detector motor position
            distance(Optional(num)): sample-detector distance
            
        Returns:
            num or dict: distance or kwargs for this function to get the distance
        """
        if detectorposition is None:
            if distance is not None:
                detectorposition = self._distancecalc(distance=distance,zerodistance=self.zerodistance)
            return {"detectorposition":detectorposition}
        else:
            return self._distancecalc(detectorposition=detectorposition,zerodistance=self.zerodistance)
            
    def calibrate_manually(self,distance):
        """Calibrate geometry based on one (distance,position) pair
        """
        distance = units.magnitude(distance,"cm")
        self.zerodistance = self._distancecalc(distance=distance,\
                            detectorposition=self.geometry.detectorposition.to("cm").magnitude)

    def calibrate_fit(self,calibrc=None,solidanglecalib=None,\
                fit=True,fixedactivearea=True,plot=False,\
                xlabel="Motor position",ylabel="Normalized Fluorescence"):
        """Calibrate geometry based in intensity vs. linear motor position (i.e. derive active area and zerodistance)
        """

        positionunits = calibrc["positionunits"]
        detectorpositionorg = units.Quantity(np.asarray(calibrc["detectorposition"]),positionunits)
        detectorposition = detectorpositionorg
        signal = np.asarray(calibrc["signal"])
        varsignal = np.asarray(calibrc["var"])
        
        dmagnitude = lambda x: x.to("cm").magnitude
        amagnitude = lambda x: x.to("cm^2").magnitude
        dquant = lambda x: units.Quantity(x,"cm")
        aquant = lambda x: units.Quantity(x,"cm^2")
        
        if solidanglecalib is None:
            # Initial parameter values
            activearea = noisepropagation.E(self.geometry.activearea)
            zerodistance = noisepropagation.E(self.zerodistance)

            distance = noisepropagation.E(self(detectorposition=detectorposition))
            rate = np.mean(signal/self.geometry.solidangle_calc(activearea=activearea,distance=distance))

            if fit:
                if fixedactivearea:
                    cactivearea = silxfit.CFIXED
                else:
                    cactivearea = silxfit.CFREE
                p0 = [rate,dmagnitude(zerodistance),amagnitude(activearea)]
                constraints = [[silxfit.CFREE,0,0],[silxfit.CFREE,0,0],[cactivearea,0,0]]

                # Fit function
                def fitfunc(x,rate,zerodistance,activearea):
                    distance = self._distancecalc(detectorposition=dquant(x),zerodistance=dquant(zerodistance))
                    sa = self.geometry.solidangle_calc(activearea=aquant(activearea),distance=distance)
                    return rate*sa
        else:
            # Fixed relationship between active area and zero-distance
            detectorpositioncalib = self.geometry.detectorposition
            def activeareafunc(zerodistance):
                distance=self._distancecalc(detectorposition=detectorpositioncalib,zerodistance=dquant(zerodistance))
                return self.geometry.solidangle_calc(distance=distance,solidangle=solidanglecalib)
            
            # Initial parameter values
            if fixedactivearea:
                czerodistance = silxfit.CFIXED
                activearea = noisepropagation.E(self.geometry.activearea)
                distance = self.geometry.solidangle_calc(activearea=activearea,solidangle=solidanglecalib)
                zerodistance = self._distancecalc(detectorposition=detectorpositioncalib,distance=distance)
            else:
                czerodistance = silxfit.CFREE
                zerodistance = noisepropagation.E(self.zerodistance)
                activearea = activeareafunc(zerodistance)
                distance = noisepropagation.E(self(detectorposition=detectorposition))
            rate = np.mean(signal/self.geometry.solidangle_calc(activearea=activearea,distance=distance))
            
            if fit:
                p0 = [rate,dmagnitude(zerodistance)]
                constraints = [[silxfit.CFREE,0,0],[czerodistance,0,0]]
                
                # Fit function
                def fitfunc(x,rate,zerodistance):
                    distance = self._distancecalc(detectorposition=dquant(x),zerodistance=dquant(zerodistance))
                    sa = self.geometry.solidangle_calc(activearea=activeareafunc(zerodistance),distance=distance)
                    return rate*sa

        if fit:
            # Fit
            p, cov_matrix, info = silxfit.leastsq(fitfunc, dmagnitude(detectorposition), signal, p0,\
                                              constraints=constraints, sigma=np.sqrt(varsignal), full_output=True)
            
            # Error and correlations:
            # H ~= J^T.J (hessian is approximated using the jacobian)
            # cov = (H)^(-1) ~= hessian
            # S = sqrt(diag(cov))
            # cor = S^(-1).cov.S^(-1)
            errors = np.sqrt(np.diag(cov_matrix))
            for i,con in enumerate(constraints):
                if con[0]==silxfit.CFIXED:
                    errors[i] = 0
            with np.errstate(divide='ignore'):
                S = np.diag(1./errors)
            cor_matrix = S.dot(cov_matrix).dot(S)
            
            # Save result
            if solidanglecalib is None:
                rate,zerodistance,activearea = p
                rate = noisepropagation.randomvariable(p[0],errors[0])
                zerodistance = noisepropagation.randomvariable(p[1],errors[1])
                activearea = noisepropagation.randomvariable(p[2],errors[2])
            else:
                rate = noisepropagation.randomvariable(p[0],errors[0])
                zerodistance = noisepropagation.randomvariable(p[1],errors[1])
                activearea = activeareafunc(zerodistance)
                
            self.zerodistance = dquant(zerodistance)
            self.geometry.detector.activearea = aquant(activearea)

        # Plot
        if plot:
            if fit:
                label = 'fit'
            else:
                label = 'current'

            plt.plot(detectorpositionorg,signal,'x',label='data')
            
            activearea = noisepropagation.E(self.geometry.activearea)
            distance = noisepropagation.E(self(detectorposition=detectorposition))
            sa = self.geometry.solidangle_calc(activearea=activearea,distance=distance)
            plt.plot(detectorpositionorg,noisepropagation.E(rate)*sa,label=label)
            
            txt = []
            if fit:
                txt.append("$\chi^2_{{red}}$ = {}".format(info["reduced_chisq"]))
            txt.append("$c$ = {:f}".format(rate))
            txt.append("$x_0$ = {:~}".format(self.zerodistance_rv.to("mm")))
            txt.append("$A$ = {:~}".format(self.geometry.activearea_rv.to("mm^2")))
            if fit:
                if solidanglecalib is None:
                    txt.append("R($c$,$x_0$) = {}".format(cor_matrix[0,1]))
                    txt.append("R($c$,$A$) = {}".format(cor_matrix[0,2]))
                    txt.append("R($x_0$,$A$) = {}".format(cor_matrix[1,2]))
                else:
                    txt.append("R($c$,$x_0$) = {}".format(cor_matrix[0,1]))
                
            off = 0.7
            for s in txt:
                plt.figtext(0.6,off,s)
                off -= 0.05
            
            ax = plt.gcf().gca()
            ax.set_xlabel("{} ({})".format(xlabel,positionunits))
            ax.set_ylabel(ylabel)
            plt.legend(loc='best')
        
    def calibrate_fit_testcorrelation2(self,signal,rate=None,zerodistance=None,detectorposition=None):

        # Fit function
        def func(x,rate,zerodistance,activearea):
            distance = x+zerodistance
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
                img[i,j] = np.sum((obs-signal)**2/signal)

        cax = plt.imshow(img, origin='lower', cmap=plt.cm.jet, interpolation='none', extent=[vrate[0],vrate[-1],vzerodistance[0],vzerodistance[-1]])
        ax = plt.gcf().gca()
        ax.set_ylabel("$x_0$ (cm)")
        ax.set_xlabel("$c_x$ (sr$^{-1}$)")
        ax.set_aspect(abs(vrate[-1]-vrate[0])/abs(vzerodistance[-1]-vzerodistance[0]))
        ax.axvline(x=rate)
        ax.axhline(y=zerodistance)
        cbar = plt.colorbar(cax,label="$\chi^2$")
        ax = plt.gcf().gca()
        
    def calibrate_fit_testcorrelation(self,signal,rate=None,zerodistance=None,detectorposition=None):

        # Fit function
        def func(x,rate,zerodistance,activearea):
            distance = x+zerodistance
            sa = self.geometry.solidangle_calc(activearea=activearea,distance=distance)
            return rate*sa

        constraints = [[silxfit.CFIXED,0,0],[silxfit.CFREE,0,0],[silxfit.CFIXED,0,0]]
        
        n = 50
        y = np.zeros(n)
        y2 = np.zeros(n)
        x = np.linspace(self.geometry.activearea-0.01,self.geometry.activearea+0.01,n)
        for i,activearea in enumerate(x):
            p0 = [rate,self.zerodistance+np.random.uniform(-1,1),activearea]
            p, cov_matrix, info = silxfit.leastsq(func, detectorposition, signal, p0,\
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


class XRFGeometry(with_metaclass(base.Centric)):

    def __init__(self,distancefunc=None,distanceargs=None,**kwargs):
        """
        Args:
            distancefunc(Optional(callable)): distance = distancefunc(**distanceargs)
            distanceargs(Optional(dict)): 
        """

        self.set_distancefunc(distancefunc)
        if distanceargs is None:
            distanceargs = {}
        self.distanceargs = distanceargs
        self._calibrcfile = None
        
        super(XRFGeometry,self).__init__(**kwargs)

    def set_distancefunc(self,distancefunc):
        if distancefunc is not None:
            distancefunc.geometry = self
        self.distancefunc = distancefunc

    @property
    def distance_rv(self):
        """Sample-detector distance in cm
        """
        if self.distancefunc is None:
            return super(XRFGeometry,self).distance
        else:
            return self.distancefunc(**self.distanceargs)
    
    @property
    def distance(self):
        return noisepropagation.E(self.distance_rv)
    
    @distance.setter
    def distance(self,distance):
        if self.distancefunc is None:
            super(XRFGeometry,self).distance = distance
        else:
            if distance is not None:
                distance = units.Quantity(distance,"cm")
            self.distanceargs.update(self.distancefunc(distance=distance))
    
    @property
    def calibrcfile(self):
        filename = self._calibrcfile
        if filename is None:
            filename = resource_filename('geometry/distancecalib_{}_{}.json'.format(self.__class__.__name__,self.detector.__class__.__name__))
        return filename
    
    @calibrcfile.setter
    def calibrcfile(self,value):
        self._calibrcfile = value
        
    def get_calibrc(self):
        filename = self.calibrcfile
        with open(filename,'r') as f:
            calibdata = json.load(f)
        return calibdata
    
    def set_calibrc(self,calibrc):
        filename = self.calibrcfile
        with open(filename,'w') as f:
            json.dump(calibrc,f,indent=2)
    
    def calibrate(self,**kwargs):
        if self.distancefunc is not None:
            if "calibrc" not in kwargs:
                kwargs["calibrc"] = self.get_calibrc()
            return self.distancefunc.calibrate_fit(**kwargs)

    def calibrate_manually(self,distance):
        if self.distancefunc is not None:
            self.distancefunc.calibrate_manually(distance)


class LinearXRFGeometry(XRFGeometry):

    def __init__(self,detectorposition=None,zerodistance=None,positionunits=None,positionsign=1,**kwargs):
        if zerodistance is not None:
            zerodistance = units.Quantity(zerodistance,positionunits)
        distancefunc = LinearMotor(zerodistance=zerodistance,positionsign=positionsign)
        super(LinearXRFGeometry,self).__init__(distancefunc=distancefunc,**kwargs)
        self.positionunits = positionunits
        self.detectorposition = detectorposition
     
    @property
    def detectorposition(self):
        return self.distanceargs["detectorposition"]
        
    @detectorposition.setter
    def detectorposition(self,value):
        if value is None:
            self.distanceargs["detectorposition"] = value
        else:
            self.distanceargs["detectorposition"] = units.Quantity(value,self.positionunits)

    @property
    def zerodistance(self):
        return self.distancefunc.zerodistance
    
    @zerodistance.setter
    def zerodistance(self,value):
        self.distancefunc.zerodistance = value
    
    def __str__(self):
        return "{}\n Detector position = {:~}".format(super(LinearXRFGeometry,self).__str__(),self.detectorposition)

    def set_calibrc(self,calibrc):
        calibrc["positionunits"] = calibrc.get("positionunits",self.positionunits)
        super(LinearXRFGeometry,self).set_calibrc(calibrc)
    
    
class sxm120(LinearXRFGeometry):

    def __init__(self,**kwargs):
        # blc10516 (April 2017)
        # detector position in mm

        kwargs["anglein"] = kwargs.get("anglein",62)
        kwargs["angleout"] = kwargs.get("angleout",49)
        kwargs["azimuth"] = kwargs.get("azimuth",0)

        kwargs["positionunits"] = kwargs.get("positionunits","mm")
        kwargs["positionsign"] = kwargs.get("sign",1)
        kwargs["zerodistance"] = units.Quantity(kwargs.get("zerodistance",56.5),"mm")
        kwargs["detectorposition"] = kwargs.get("detectorposition",0)

        super(sxm120,self).__init__(**kwargs)
        

class sxm90(LinearXRFGeometry):

    def __init__(self,**kwargs):
        # blc10516 (April 2017)
        # detector position in mm
        
        kwargs["anglein"] = kwargs.get("anglein",62)
        kwargs["angleout"] = kwargs.get("angleout",28)
        kwargs["azimuth"] = kwargs.get("azimuth",0)

        kwargs["positionunits"] = kwargs.get("positionunits","mm")
        kwargs["positionsign"] = kwargs.get("sign",1)
        kwargs["zerodistance"] = units.Quantity(kwargs.get("zerodistance",85.571),"mm")
        kwargs["detectorposition"] = kwargs.get("detectorposition",0)
        
        super(sxm90,self).__init__(**kwargs)


class microdiff(LinearXRFGeometry):

    def __init__(self,**kwargs):
        # blc10516 (April 2017)
        # detector position in mm
        
        kwargs["anglein"] = kwargs.get("anglein",62)
        kwargs["angleout"] = kwargs.get("angleout",28)
        kwargs["azimuth"] = kwargs.get("azimuth",0)

        kwargs["positionunits"] = kwargs.get("positionunits","mm")
        kwargs["positionsign"] = kwargs.get("sign",-1)
        kwargs["zerodistance"] = units.Quantity(kwargs.get("zerodistance",-56.667),"mm")
        kwargs["detectorposition"] = kwargs.get("detectorposition",10.)
        
        super(microdiff,self).__init__(**kwargs)
                        

factory = XRFGeometry.factory

