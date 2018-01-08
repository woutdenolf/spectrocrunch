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
import numbers

from .. import xraylib
from .. import ureg
from ..common import instance
from ..math import fit1d
from ..common import hashable
from ..common import listtools
from ..common.Enum import Enum

class FluoLine(hashable.Hashable):

    @staticmethod
    def getlinename(line):
        # Return IUPAC instead of SIEGBAHN when possible
        if isinstance(line,FluoLine):
            return line.name
        elif isinstance(line,numbers.Number):
            names = xraylib._code_to_line[line]
            
            n = len(names)
            if n==1:
                return names[0]
            elif n==2:
                n0 = len(names[0])
                n1 = len(names[1])
                if n0==n1:
                    candidate = names[0]
                    if names[0][1]=='A' or names[0][1]=='B':
                        return names[1]
                    else:
                        return names[0]
                elif n0>n1:
                    return names[0]
                else:
                    return names[1]
            else:
                return names[np.argmax([len(name) for name in names])]
        else:
            if line not in xraylib._line_to_code.keys():
                raise ValueError("Unknown line name {}".format(line))
            return line

    @staticmethod
    def getlinecode(line):
        if isinstance(line,FluoLine):
            return line.code
        elif isinstance(line,numbers.Number):
            if line not in xraylib._code_to_line.keys():
                raise ValueError("Unknown line code {}".format(line))
            return line
        else:
            return xraylib._line_to_code[line]
    
    @staticmethod
    def decompose(code):
        """
        Args:
            code(int):
        Return:
            list(int): sub lines
        """
        if code in xraylib._composites:
            return xraylib._composites[code]
        else:
            return [code]
        
    @classmethod
    def all_lines(cls):
        return list(set(cls(code) for line in range(xraylib._line_max,xraylib._line_min-1,-1) for code in cls.decompose(line)))
    
    @classmethod
    def factory(cls,shells=None,fluolines=None,energybounds=None):
        """Generate FluoLine instances, possibly limited by shell or energy range

        Args:
            shells(Optional(array(int or str or Shell))): 
            fluolines(Optional(array(int or str or FluoLine))):  None or [] -> explicite all
            energybounds(Optional(3-tuple)): element or num or str, lower energy bound, higher energy bound
            
        Returns:
            list(FluoLine)
        """

        # All lines or the lines given
        if fluolines is None:
            lines = cls.all_lines()
        else:
            if instance.isarray(fluolines):
                if fluolines==[]:
                    lines = cls.all_lines()
                else:
                    lines = [FluoLine(line) for line in fluolines]
            else:
                lines = [FluoLine(fluolines)]
                
        # Shell selection
        if shells is not None:
            if instance.isarray(shells):
                shellnames = [Shell.getshellname(s) for s in shells]
            else:
                shellnames = [Shell.getshellname(shells)]
            valid = lambda linecode: any(any(cls.getlinename(code).startswith(shellname) for shellname in shellnames) for code in cls.decompose(linecode)) 
            lines = [line for line in lines if valid(line.code)]

        # Energy selection
        if energybounds is not None:
            Z = energybounds[0]
            valid = lambda energy: energy >= energybounds[1] and energy <= energybounds[2] and energy!=0
            lines = [line for line in lines if valid(line.energy(Z))]
                
        # Return
        return lines
            
    def __init__(self,line):
        """
        Args:
            code(int or str): xraylib line
        """
        if isinstance(line,self.__class__):
            self.code = line.code
            self.name = line.name
        else:
            self.code = self.getlinecode(line)
            self.name = self.getlinename(line)

    def _cmpkey(self):
        """For comparing and sorting
        """
        return self.code

    def _stringrepr(self):
        """Unique representation of an instance
        """
        return self.name

    def radrate(self,Z):
        """Radiative rate of a line: probability of this line / probabilty of fluorescence
        
        Args:
            Z(num): atomic number
            
        Returns:
            num
        """
        rate = xraylib.RadRate(Z,self.code)
        
        if rate==0:
            # Maybe the radiative rate of one of its composites is known
            if self.code in xraylib._rcomposites:
                for comp in xraylib._rcomposites[self.code]:
                    if comp>=0:
                        continue
                    rate = xraylib.RadRate(Z,comp)
                    if rate!=0:
                        # In abscence of a better assumption, assume all lines
                        # in the composite are equally probable
                        return rate/len(xraylib._composites[comp])
                        
        return rate
    
    def energy(self,Z):
        """
        Args:
            Z(num): atomic number
        Returns:
            num
        """
        energy = xraylib.LineEnergy(Z,self.code)
        
        if energy==0:
            # Maybe the energy of one of its composites is known
            if self.code in xraylib._rcomposites:
                for comp in xraylib._rcomposites[self.code]:
                    if comp>=0:
                        continue
                    energy = xraylib.LineEnergy(Z,comp)
                    if energy!=0:
                        return energy
                        
        return energy

    @property
    def shell(self):
        """Shell to which this line belongs to
        """
        for shellname in xraylib._shell_to_code:
            if self.name.startswith(shellname):
                return Shell(shellname)
        raise RuntimeError("Cannot find the shell of fluorescence line {}".format(self))


class Shell(hashable.Hashable):

    @staticmethod
    def getshellname(shell):
        if isinstance(shell,Shell):
            return shell.name
        elif isinstance(shell,numbers.Number):
            return xraylib._code_to_shell[shell]
        else:
            if shell not in xraylib._shell_to_code.keys():
                raise ValueError("Unknown shell name {}".format(shell))
            return shell
            
    @staticmethod
    def getshellcode(shell):
        if isinstance(shell,Shell):
            return shell.code
        elif isinstance(shell,numbers.Number):
            if shell not in xraylib._code_to_shell.keys():
                raise ValueError("Unknown shell code {}".format(shell))
            return shell
        else:
            return xraylib._shell_to_code[shell]

    @classmethod
    def all_shells(cls,fluolines=None):
        shells = range(xraylib._shell_min,xraylib._shell_max+1)
        return [cls(shell,fluolines=fluolines) for shell in shells]       

    @classmethod
    def factory(cls,energybounds=None):
        if energybounds is None:
            return cls.all_shells()
        else:
            Z = energybounds[0]
            valid = lambda energy: energy >= energybounds[1] and energy <= energybounds[2]
            shells = cls.all_shells()
            shells = [s for s in shells if valid(s.edgeenergy(Z))]
            for s in shells:
                s._fluolines = FluoLine.factory(shells=[s],energybounds=energybounds)
            shells = [s for s in shells if s._fluolines]
            return shells
            
    @classmethod
    def pymcafactory(cls,energybounds=None):
        return list(set("".join(re.split("[^a-zA-Z]*", str(s))) for s in cls.factory(energybounds=energybounds)))

    def __init__(self,shell,fluolines=None):
        """
        Args:
            code(int or str or Shell): shell code or name
            fluolines(Optional(array(int or str or Fluoline))): emission lines
        """
        self.code = self.getshellcode(shell)
        self.name = self.getshellname(shell)
        if isinstance(shell,Shell) and fluolines is None:
            self._fluolines = shell._fluolines
        else:
            if fluolines is None:
                self._fluolines = None # all lines
            else:
                self._fluolines = FluoLine.factory(shells=[self],fluolines=fluolines)

    @property
    def fluolines(self):
        if self._fluolines is None:
            # None means all lines, but I need them explicitely now
            self._fluolines = FluoLine.factory(shells=[self])
        return self._fluolines

    def _cmpkey(self):
        """For comparing and sorting
        """
        return self.code

    def _stringrepr(self):
        """Unique representation of an instance
        """
        return self.name

    def fluoyield(self,Z):
        """Fluorescence yield for this shell: probability for fluorescence / probability of shell ionization
        """
        return xraylib.FluorYield(Z,self.code)

    def radrate(self,Z):
        """Radiative rate of a shell: probabilities of selected lines / probabilty of fluorescence 
        """
        if self._fluolines is None:
            return [1] # ALL lines
        else:
            return [l.radrate(Z) for l in self._fluolines]

    def partial_fluoyield(self,Z,decomposed=False):
        """Probability of selected lines / probability of shell ionization
        """
        if decomposed:
            fluoyield = self.fluoyield(Z)
            return {l:fluoyield*l.radrate(Z) for l in self.fluolines}
        else:
            if self._fluolines is None:
                return self.fluoyield(Z)
            else:
                return self.fluoyield(Z)*sum(self.radrate(Z))

    def edgeenergy(self,Z):
        return xraylib.EdgeEnergy(Z,self.code)


class FluoZLine(hashable.Hashable):

    def __init__(self,element,line):
        self.line = line
        self.element = element
        
    def __getattr__(self,attr):
        try:
            return getattr(self.line,attr)
        except:
            return getattr(self.element,attr)

    def _cmpkey(self):
        """For comparing and sorting
        """
        return self._stringrepr()

    def _stringrepr(self):
        """Unique representation of an instance
        """
        return "{}-{}".format(self.element,self.line)
        
    def energy(self,**kwargs):
        return self.line.energy(self.element.Z)


class RayleighLine(hashable.Hashable):

    def __init__(self,energysource):
        self.energysource = energysource
        self.name = "Rayleigh"
        
    def _cmpkey(self):
        """For comparing and sorting
        """
        return self.name

    def _stringrepr(self):
        """Unique representation of an instance
        """
        return self.name

    def energy(self,**kwargs):
        return self.energysource
        
        
class ComptonLine(hashable.Hashable):

    def __init__(self,energysource):
        self.energysource = energysource
        self.name = "Compton"

    def _cmpkey(self):
        """For comparing and sorting
        """
        return self.name

    def _stringrepr(self):
        """Unique representation of an instance
        """
        return self.name

    def energy(self,scatteringangle=0):
        if scatteringangle==0:
            return self.energysource
        delta = ureg.Quantity(1-np.cos(np.radians(scatteringangle)),"1/(m_e*c^2)").to("1/keV","spectroscopy").magnitude
        return self.energysource/(1+self.energysource*delta) 
        
        
class Spectrum(object):

    TYPES = Enum(['crosssection','interactionyield'])
    # crosssection: cm^2/g
    # interactionyield: crosssection.density.thickness (dimensionless)

    def __init__(self):
        self.cs = {}
        self.density = None
        self.xlim = None
        self.title = None
        self.xlabel = "Energy (keV)"
        self.type = None
   
    @property
    def lines(self):
        for line,cs in self.cs.items():
            if np.sum(cs)==0:
                continue
            yield line

    @property
    def sortedcs(self):
        for line,cs in sorted(self.cs.items(),key=lambda x: max(instance.asarray(x[0].energy(scatteringangle=0.1)))):
            if np.sum(cs)==0:
                continue
            yield line,cs

    @property
    def probabilities(self):
        """Interaction probability per cm and per srad
        """
        for line,cs in self.cs.items():
            if np.sum(cs)==0:
                continue
            yield line,cs*(self.density/4*np.pi) #TODO: assume all isotropic for now

    def energies(self,**kwargs):
        for line in self.cs:
            energy = instance.asarray(line.energy(**kwargs))
            for en in energy:
                yield line,en

    def expandedenergies(self,**kwargs):
        for line,cs in self.cs:
            try:
                n = len(cs)
            except:
                n = 1
            energy = instance.asarray(line.energy(**kwargs))
            if n>1 and len(energy)==1:
                energy = np.repeat(energy,n)
            return line,energy

    def plotlineprobability(self,energy,probability,detector=None,**kwargs):
        if detector is None:
            plt.plot([energy,energy],[0,probability],**kwargs)
            h = probability
        else:
            FWHM = detector.linewidth(energy)
            sx = FWHM/(2*np.sqrt(2*np.log(2)))
            k = 4
            x = np.linspace(energy-k*sx,energy+k*sx,50)
            y = fit1d.gaussian(x,energy,sx,probability)
            h = max(y)
            plt.plot(x,y,**kwargs)
        return h

    def plotprobability(self,detector=None,out=False,mark=True,log=False):
        ax = plt.gca()
        
        colors = {}
        markers = {}
        bmark = False
        
        if detector is None:
            scatteringangle = 90.
        else:
            scatteringangle = detector.scatteringangle

        for line,cs in self.sortedcs:
                
            element = getattr(line, 'element', None)
            energy = instance.asarray(line.energy(scatteringangle=scatteringangle))
            probability = cs*self.density

            # Scattering or fluorescence
            if element is None: # scattering
                key = str(line)
            
                # Line label + color
                label = [None]*len(energy)
                label[-1] = key
                color = next(ax._get_lines.prop_cycler)['color']
                colors[key] = color
                
                # Line marker
                if mark:
                    markers[key] = {"name":key,"energy":energy,"probability":probability,"height":probability}
                    bmark = True
                    
                if log:
                    print(line,energy,probability)
            else: # fluorescence
                key = "{} {}".format(element,line.shell)

                # Line label + color
                label = [None]*len(energy)
                if key in colors:
                    color = colors[key]
                else:
                    label[-1] = key
                    color = next(ax._get_lines.prop_cycler)['color']
                    colors[key] = color

                # Key marker
                if mark:
                    if key in markers:
                        bmark = sum(probability)>sum(markers[key]["probability"])
                    else:
                        bmark = True
                    if bmark:
                        markers[key] = {"name":str(line),"energy":energy,"probability":probability,"height":probability}

                if log:
                    print(line,energy,probability)
            
            height = [self.plotlineprobability(en,prob,color=color,label=lab,detector=detector) for en,prob,lab in zip(energy,probability,label)]
            
            if bmark:
                ind = np.argmax(height)
                markers[key]["energy"] = markers[key]["energy"][ind]
                markers[key]["height"] = height[ind]

        for key,marker in markers.items():
            plt.annotate(marker["name"], xy=(marker["energy"], marker["height"]),color=colors[key])

        if detector is not None and log:
            ax.set_yscale('log', basey=10)
            plt.ylim([1,None])
            
        plt.legend(loc='best')
        plt.xlabel(self.xlabel)
        if self.type==self.TYPES.crosssection:
            plt.ylabel("Probability (1/cm)")
        else:
            plt.ylabel("Generation yield")
        plt.xlim(self.xlim)
        try:
            plt.title(self.title)
        except UnicodeDecodeError:
            pass


