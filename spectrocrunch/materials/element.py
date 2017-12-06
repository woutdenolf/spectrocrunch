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

import numpy as np
import os
import tempfile
import time
from scipy import interpolate
import json
import numbers
import re
try:
    import iotbx.cif
except:
    pass

import fdmnes

from ..common.hashable import Hashable
from ..common.instance import isarray
from .. import xraylib
from .. import ureg

def parseelementsymbol(symb):
    if isinstance(symb,str):
        if symb.isdigit():
            Z = int(symb)
            name = xraylib.AtomicNumberToSymbol(Z)
        else:
            Z = xraylib.SymbolToAtomicNumber(symb.title())
            name = symb
    elif isinstance(symb,numbers.Integral):
        Z = symb
        name = xraylib.AtomicNumberToSymbol(Z)
    elif isinstance(symb,Element):
        Z,name = symb.Z,symb.name
    else:
        raise ValueError("Unknown element symbol or number")
    return Z,name  

class FluoLine(Hashable):

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
            fluolines(Optional(array(int or str or FluoLine))):  None or [] -> explicit all
            energybounds(Optional(3-tuple)): element or num or str, lower energy bound, higher energy bound
            
        Returns:
            list(FluoLine)
        """

        # All lines or the lines given
        if fluolines is None:
            lines = cls.all_lines()
        else:
            if isarray(fluolines):
                if fluolines==[]:
                    lines = cls.all_lines()
                else:
                    lines = [FluoLine(line) for line in fluolines]
            else:
                lines = [FluoLine(fluolines)]
                
        # Shell selection
        if shells is not None:
            if isarray(shells):
                shellnames = [Shell.getshellname(s) for s in shells]
            else:
                shellnames = [Shell.getshellname(shells)]
            valid = lambda linecode: any(any(cls.getlinename(code).startswith(shellname) for shellname in shellnames) for code in cls.decompose(linecode)) 
            lines = [line for line in lines if valid(line.code)]

        # Energy selection
        if energybounds is not None:
            Z,_ = parseelementsymbol(energybounds[0])
            valid = lambda energy: energy >= energybounds[1] and energy <= energybounds[2]
            lines = [line for line in lines if valid(line.lineenergy(Z))]
        
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
    
    def lineenergy(self,Z):
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

    def shell(self):
        """Shell to which this line belongs to
        """
        for shellname in xraylib._shell_to_code:
            if self.name.startswith(shellname):
                return Shell(shellname)
        raise RuntimeError("Cannot find the shell of fluorescence line {}".format(self))


class Shell(Hashable):

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
        alls = cls.all_shells()
        if energybounds==0:
            return alls
        else:
            Z,_ = parseelementsymbol(energybounds[0])
            valid = lambda energy: energy >= energybounds[1] and energy <= energybounds[2]
            return [s for s in alls if valid(s.edgeenergy(Z))]
    
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
            return FluoLine.factory(shells=[self])
        else:
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

    def partial_fluoyield(self,Z):
        """Probability of selected lines / probability of shell ionization
        """
        if self._fluolines is None:
            return self.fluoyield(Z)
        else:
            return self.fluoyield(Z)*sum(self.radrate(Z))

    def edgeenergy(self,Z):
        return xraylib.EdgeEnergy(Z,self.code)
        
    
class Element(Hashable):
    """Interface to chemical elements
    """

    def __init__(self,symb):
        """
        Args:
            symb(int or str or Element)
        """
        
        # Atomic number
        self.Z,self.name = parseelementsymbol(symb)
        
        # Atomic properties
        self.MM = xraylib.AtomicWeight(self.Z)

        # Information for calculating partial cross-sections
        self.shells = []
        
    def _cmpkey(self):
        """For comparing and sorting
        """
        return -self.Z

    def _stringrepr(self):
        """Unique representation of an instance
        """
        return self.name

    def markabsorber(self,symb,shells=None,fluolines=None):
        """Marking an element's shells and lines has following effect:
            - partial absorption cross-section is not zero
            - when no fluolines are given: all lines for each shell are taken into account

        Args:
            symb(str): element symbol
            shells(Optional(array(int))): None -> all
            fluolines(Optional(array(int))): None -> all, [] -> explicit all
        """
        if self.name==symb:
            # Shell names
            if shells is None:
                self.shells = Shell.all_shells(fluolines=fluolines)
            else:
                if isarray(shells):
                    self.shells = [Shell(shell,fluolines=fluolines) for shell in shells]
                else:
                    self.shells = [Shell(shells,fluolines=fluolines)]

    def unmarkabsorber(self):
        self.shells = []

    def isabsorber(self):
        return len(self.shells)!=0

    @property
    def density(self):
        return xraylib.ElementDensity(self.Z)

    def molarmass(self):
        return self.MM
        
    def weightfractions(self):
        return dict([(self,1.)])

    def molefractions(self,total=True):
        return dict([(self,1.)])
    
    @property
    def nelements(self):
        return 1

    def _get_fdmnes_energyrange(self,Eabs,edgeenergy,decimals=6):
        """Calculate energies based on boundaries and step sizes:
            energyrange = [E0,step0,E1,step1,E2] (eV, relative to the edge)
           Eabs in keV
        """

        E = (Eabs - edgeenergy)*1000 # absolute (keV) to relative (eV)

        dE = np.around(E[1:]-E[0:-1],decimals=decimals-3)
        ind = np.argwhere(dE[1:]-dE[0:-1])
        nblocks = len(ind)+1

        energyrange=np.empty(2*nblocks+1,dtype=E.dtype)
        energyrange[0] = E[0]
        energyrange[-1] = E[-1]

        if nblocks==1:
            energyrange[1] = dE[0]
        else:
            inddest = np.arange(1,2*nblocks-2,2)
            energyrange[inddest] = dE[ind]
            energyrange[inddest+1] = E[ind+1]
            energyrange[-2] = dE[ind[-1]+1]

        # enlarge the fdmnes range a bit
        add = 5 # eV
        energyrange[0] -= np.ceil(add/energyrange[1])*energyrange[1]
        energyrange[-1] += np.ceil(add/energyrange[-2])*energyrange[-2]

        return energyrange

    def mass_att_coeff(self,E,environ=None,decimals=6,refresh=False):
        """Mass attenuation coefficient (cm^2/g, E in keV). Use for transmission XAS.

        Args:
            E(num or array-like): energy (keV)
            environ(dict): chemical environment of this element
            decimals(Optional(num)): precision of energy in keV
            refresh(Optional(bool)): force re-simulation if used

        Returns:
            num or np.array
        """
        bnum = not isarray(E)
        if bnum:
            E = [E]
        cs = np.empty(len(E),dtype=np.float64)

        # Total
        for i in range(len(E)):
            cs[i] = xraylib.CS_Total_Kissel(self.Z,np.float64(E[i]))

        # Replace part by simulation
        cs = self._replace_partial_mass_abs_coeff(cs,E,environ=environ,decimals=decimals,refresh=refresh)

        if bnum:
            return cs[0]
        else:
            return cs

    def mass_abs_coeff(self,E,environ=None,decimals=6,refresh=False):
        """Mass absorption coefficient (cm^2/g, E in keV).

        Args:
            E(num or array-like): energy (keV)
            environ(dict): chemical environment of this element
            decimals(Optional(num)): precision of energy in keV
            refresh(Optional(bool)): force re-simulation if used

        Returns:
            num or np.array
        """
        bnum = not isarray(E)
        if bnum:
            E = [E]
        cs = np.empty(len(E),dtype=np.float64)

        # Total
        for i in range(len(E)):
            cs[i] = xraylib.CS_Photo_Total(self.Z,np.float64(E[i]))

        # Replace part by simulation
        cs = self._replace_partial_mass_abs_coeff(cs,E,environ=environ,decimals=decimals,refresh=refresh)

        if bnum:
            return cs[0]
        else:
            return cs

    def _replace_partial_mass_abs_coeff(self,cs,E,environ=None,decimals=6,refresh=False):
        """
        """
        if environ is not None:
            # Subtract partial cross-sections (summed over selected shells)
            cs -= sum(self._CS_Photo_Partial_DB(E).values())

            ind = np.argwhere(cs<0)
            if len(ind)>0:
                ind2 = np.argwhere(cs>=0)
                f = interpolate.interp1d(E[ind2].flatten(),cs[ind2].flatten(),bounds_error=False)
                cs[ind] = f(E[ind])

            # Add partial cross-sections (summed over selected shells)
            cs += sum(self._CS_Photo_Partial_SIM(E,environ,decimals=decimals,refresh=refresh).values())

        return cs

    def partial_mass_abs_coeff(self,E,environ=None,decimals=6,refresh=False):
        """Mass absorption coefficient for the selected shells (cm^2/g, E in keV).

        Args:
            E(num or array-like): energy (keV)
            environ(dict): chemical environment of this element
            decimals(Optional(num)): precision of energy in keV
            refresh(Optional(bool)): force re-simulation if used

        Returns:
            num or array: sum_S[tau(E,S)]
        """
        bnum = not isarray(E)
        if bnum:
            E = [E]

        if not self.isabsorber():
            if bnum:
                return np.float64(0)
            else:
                return np.zeros(len(E),dtype=np.float64)

        if environ is None:
            cs = self._CS_Photo_Partial_DB(E)
        else:
            cs = self._CS_Photo_Partial_SIM(E,environ,decimals=decimals,refresh=refresh)

        # sum over selected shells
        cs = sum(cs.values())
        if bnum:
            return cs[0]
        else:
            return cs

    def fluorescence_cross_section(self,E,environ=None,decimals=6,refresh=False,decomposed=False):
        """XRF cross section for the selected shells and lines (cm^2/g, E in keV). Use for fluorescence XAS.

            muXRF(E) = sum_{S}[tau(E,S)*fluoyield(S)*sum_{L}[radrate(S,L)]]

        Args:
            E(num or array-like): energy (keV)
            environ(dict): chemical environment of this element
            decimals(Optional(num)): precision of energy in keV
            refresh(Optional(bool)): force re-simulation if used
            decomposed(Optional(bool)): output as dictionary

        Returns:
            num or np.array: muXRF
            dict: S:tau(E,S)
        """
        bnum = not isarray(E)
        if bnum:
            E = [E]

        if not self.isabsorber():
            if bnum:
                return np.float64(0)
            else:
                return np.zeros(len(E),dtype=np.float64)

        if environ is None:
            cs = self._CS_Photo_Partial_DB(E)
        else:
            cs = self._CS_Photo_Partial_SIM(E,environ,decimals=decimals,refresh=refresh)

        if decomposed:
            return cs
        else:
            # sum over selected shells, weighted but the total fluoyield of the selected lines
            cs = sum([c*shell.partial_fluoyield(self.Z) for shell,c in cs.items()])
            if bnum:
                return cs[0]
            else:
                return cs

    def _xraylib_method(self,method,E):
        method = getattr(xraylib,method)
        if isarray(E):
            ret = np.empty(len(E),dtype=np.float64)

            for i in range(len(E)):
                ret[i] = method(self.Z,np.float64(E[i]))
        else:
            ret = method(self.Z,np.float64(E))

        return ret
        
    def scattering_cross_section(self,E,environ=None,decimals=6,refresh=False):
        """Scattering cross section (cm^2/g, E in keV).

        Args:
            E(num or array-like): energy (keV)
            environ(dict): chemical environment of this element
            decimals(Optional(num)): precision of energy in keV
            refresh(Optional(bool)): force re-simulation if used

        Returns:
            num or np.array
        """

        return self._xraylib_method("CS_Rayl",E)+self._xraylib_method("CS_Compt",E)

    def rayleigh_cross_section(self,E,environ=None,decimals=6,refresh=False):
        """Rayleigh cross section (cm^2/g, E in keV).

        Args:
            E(num or array-like): energy (keV)
            environ(dict): chemical environment of this element
            decimals(Optional(num)): precision of energy in keV
            refresh(Optional(bool)): force re-simulation if used

        Returns:
            num or np.array
        """
        return self._xraylib_method("CS_Rayl",E)

    def compton_cross_section(self,E,environ=None,decimals=6,refresh=False):
        """Rayleigh cross section (cm^2/g, E in keV).

        Args:
            E(num or array-like): energy (keV)
            environ(dict): chemical environment of this element
            decimals(Optional(num)): precision of energy in keV
            refresh(Optional(bool)): force re-simulation if used

        Returns:
            num or np.array
        """
        return self._xraylib_method("CS_Compt",E)

    def scatfact_classic_re(self,E,theta=None,environ=None,decimals=6,refresh=False):
        """Real part of atomic form factor

        Args:
            E(num or array-like): energy (keV)
            theta(Optional(num or array-like)): scattering angle (rad)
            environ(dict): chemical environment of this element
            decimals(Optional(num)): precision of energy in keV
            refresh(Optional(bool)): force re-simulation if used

        Returns:
            num or np.array
        """
        if theta is None:
            return self.Z
        else:
            #q = Q/4.pi
            q = np.sin(theta/2)/ureg.Quantity(E,'keV').to("angstrom","spectroscopy").magnitude
            return self._xraylib_method("FF_Rayl",q)
            
    def scatfact_re(self,E,theta=None,environ=None,decimals=6,refresh=False):
        """Real part of atomic form factor

        Args:
            E(num or array-like): energy (keV)
            theta(Optional(num or array-like)): scattering angle (degrees)
            environ(dict): chemical environment of this element
            decimals(Optional(num)): precision of energy in keV
            refresh(Optional(bool)): force re-simulation if used

        Returns:
            num or np.array
        """
        return self.scatfact_classic_re(E,theta=theta,environ=environ,decimals=decimals,refresh=refresh) + self._xraylib_method("Fi",E)
    
    def scatfact_im(self,E,environ=None,decimals=6,refresh=False):
        """Imaginary part of atomic form factor

        Args:
            E(num or array-like): energy (keV)
            environ(dict): chemical environment of this element
            decimals(Optional(num)): precision of energy in keV
            refresh(Optional(bool)): force re-simulation if used

        Returns:
            num or np.array
        """
        return self._xraylib_method("Fii",E)
    
    def _get_multiplicity(self,struct):
        scat = struct.scatterers()
        ret = 0.
        for s in scat:
            if s.scattering_type==self.name:
                ret += s.occupancy*s.multiplicity()
        return ret

    def _CS_Photo_Partial_DB(self,E):
        """Get the partial photoionization cross section from xraylib (E in keV).

        Args:
            E(array): energy (keV)

        Returns
            dict: Shell:tau(E,Shell)
        """
        cs = {}
        for shell in self.shells:
            cs[shell] = np.empty(len(E),dtype=np.float64)
            for i in range(len(E)):
                cs[shell][i] = xraylib.CS_Photo_Partial(self.Z,shell.code,E[i])

        return cs

    def _CS_Photo_Partial_SIM(self,E,environ,decimals=6,refresh=False,fluo=True):
        """Calculate the partial photoionization cross section with fdmnes (E in keV).

        Args:
            E(array): energy (keV)
            environ(dict): chemical environment of this element
            decimals(Optional(num)): precision of energy in keV
            refresh(Optional(bool)): force re-simulation if used

        Returns
            dict: Shell:tau(E,Shell)
        """

        # Initialize simulation
        sim = fdmnes.fdmnes(environ.ciffile, resonant=self.name)

        # Energy range
        sim.P.Energpho = False # relative energies as output

        # Simultation settings
        #sim.P.Radius = 3.5 # Radius of the cluster for calculation
        sim.P.Radius = 7.
        sim.P.Rpotmax = sim.P.Radius + 5 # Radius of the cluster for potential calculation
        sim.P.Quadrupole = True # multipole approximation
        sim.P.Green = False # MS instead of FDM (faster but less accurate)
        sim.P.TDDFT = True # multi electron correction
        sim.P.Convolution = True # 
        sim.P.Density = False # save density of states

        # Do simulation
        cs = {}

        for shell in self.shells:
            cs[shell] = np.empty(len(E),dtype=np.float64)

            # Select edge
            sim.P.Edge = shell.name
            filebase = os.path.splitext(os.path.basename(environ.ciffile))[0]+"_"+self.name+"_"+sim.P.Edge

            # Relative energy range in eV
            Range = self._get_fdmnes_energyrange(E,shell.edgeenergy(self.Z),decimals=decimals)
            sim.P.Range = tuple(Range)

            # Read previous configuration file
            fcfg = os.path.join(tempfile.gettempdir(),filebase+".json")
            if os.path.isfile(fcfg) and refresh==False:
                with open(fcfg,'r') as f:
                    config = json.load(f)
            else:
                config = {}

            # Write simulation file
            finput = os.path.join(tempfile.gettempdir(),filebase+".txt")
            sim.WriteInputFile(finput, overwrite=True)

            # Perform simulation if not already done
            if 'Range' in config:
                Range_prev = np.array(config["Range"])
            else:
                Range_prev = Range*0
            config["Range"] = sim.P.Range

            if len(Range_prev)==len(Range):
                brun = not np.allclose(Range_prev,Range)
            else:
                brun = True
            if brun:
                sim.Run(wait=False)
                while True:
                    NumRunning = sum([sim.Status(i)==False for i in range(len(sim.proc))])
                    if NumRunning == 0:
                        break
                    time.sleep(5)

            # Get data
            data = sim.get_XANES(conv = True)

            # Convert data to absorbance
            data[:,0] /= 1000. # ev -> keV
            data[:,0] += shell.edgeenergy(self.Z) # rel -> abs
            data[:,1] *= xraylib.AVOGNUM*1E6/self.MM # Absorption cross section (Mbarn) -> mass absorption coefficient (cm^2/g)

            # Element multiplicity in the unit cell
            if "nmult" in config:
                nmult = config["nmult"]
            else:
                # TODO: cctbx is too heavy (information is also in the bav file)
                nmult = self._get_multiplicity(environ.structure)
                #nmult = np.round(data[-1,1]/xraylib.CS_Photo_Partial(self.Z,Shell,E[-1])) # not precise enough
                config["nmult"] = nmult
            data[:,1] /= nmult

            # Energy interpolation and keep
            f = interpolate.interp1d(data[:,0],data[:,1],bounds_error=False)
            cs[shell] = f(E)
        
            # Write configuration file
            with open(fcfg,'w') as f:
                json.dump(config,f,indent=2)

        return cs

    def get_energy(self,energyrange,defaultinc=1):
        """Get absolute energies (keV) from a relative energy range (eV)
        """
        nshells = len(self.shells)
        if nshells==0:
            return None
        elif nshells==1:
            return self.shells[0].edgeenergy(self.Z) + self._get_fdmnes_energies(energyrange)/1000.
        else:
            return self._get_fdmnes_energies(self._get_absolute_energyrange(energyrange,defaultinc=defaultinc))

    def _get_fdmnes_energies(self,energyrange):
        """Calculate energies based on boundaries and step sizes:
            energyrange = [E0,step0,E1,step1,E2]
        """

        # Number of steps in each region
        nblocks = len(energyrange)//2
        nsteps = np.empty(nblocks,dtype=int)
        e = energyrange[0]
        for i in range(nblocks):
            b = e
            e = energyrange[2*i+2]
            inc = energyrange[2*i+1]
            nsteps[i] = np.ceil((e-b)/inc)
            e = b + nsteps[i]*inc
            if i == nblocks-1:
                if e > energyrange[-1]:
                    nsteps[i]-=1

        ret = np.empty(nsteps.sum()+1,dtype=energyrange.dtype)
        ret[0] = energyrange[0]
        off = 1
        for i in range(nblocks):
            inc = energyrange[2*i+1]
            ret[off:off+nsteps[i]] = ret[off-1] + np.arange(1,nsteps[i]+1)*inc
            off += nsteps[i]

        return ret

    def _get_absolute_energyrange(self,energyrange,defaultinc=1):
        """Convert relative energy range (eV) to absolute energy range (keV)
        """

        # Boundaries
        nblocks = len(energyrange)//2
        indE = np.arange(0,len(energyrange),2)
        nshells = len(self.shells)
        nbounds = len(indE)
        energies = np.empty(nshells*nbounds,dtype=energyrange.dtype)
        for i in range(nshells):
            energies[i*nbounds:(i+1)*nbounds] = self.shells[i].edgeenergy(self.Z) + energyrange[indE]/1000.

        # Put unique boundaries in new range array
        E = np.unique(energies)
        newnblocks = len(E)-1
        newrange = np.empty(2*newnblocks+1,dtype=energyrange.dtype)
        newindE = np.arange(0,len(newrange),2)
        newindD = np.arange(1,2*newnblocks,2)
        newrange[newindE] = E
        newrange[newindD] = defaultinc/1000.
        
        # Determine step sizes
        for j in range(newnblocks):
            m = newrange[2*j] + (newrange[2*j+2]-newrange[2*j])/2
            inc = 0
            
            for i in range(nshells):
                b = i*nbounds
                e = b+nbounds-1
                if m > energies[b] and m < energies[e]:
                    ubound = 2*np.argmax((energies[b:e+1]-m)>0)
                    d = energyrange[ubound-1]/1000.
                    if inc==0:
                        inc = d
                    else:
                        inc = min(inc,d)

            if inc != 0:
                newrange[2*j+1] = inc

        return newrange

    def pymcaformat(self):
        r = self.weightfractions()
        value = {'Comment': self.name,
                'CompoundFraction': [1.],
                'Thickness': 0.,
                'Density': self.density,
                'CompoundList': ['{}1'.format(self.name)]}
        return self.name,value
        
