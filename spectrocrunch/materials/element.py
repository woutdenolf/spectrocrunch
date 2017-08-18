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

import xraylib
xraylib.XRayInit()
xraylib.SetErrorMessages(0)

import fdmnes

import numpy as np
import os
import tempfile
import time
from scipy import interpolate
import json
try:
    import iotbx.cif
except:
    pass

from ..common.hashable import Hashable

class fluoline(Hashable):
    _all = {s.split('_')[0]:xraylib.__dict__[s] for s in xraylib.__dict__.keys() if s.endswith("_LINE")}
    
    def __init__(self,code):
        """
        Args:
            code(int): xraylib line code
        """
        self.code = code
        self.name = fluoline.getname(self.code)

    def _cmpkey(self):
        """For comparing and sorting
        """
        return self.code

    def _stringrepr(self):
        """Unique representation of an instance
        """
        return self.name

    def radrate(self,Z):
        return xraylib.RadRate(Z,self.code)

    @classmethod
    def getname(cls,code):
        """IUPAC name
        """
        names = [name for name,c in cls._all.items() if c == code]

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
    
    @classmethod
    def factory(cls,codes,shellname):
        """Generate fluoline instances from codes
            - only from given shell
            - decompose composite lines
            - remove duplicates

        Args:
            codes(array(int)): xraylib line codes

        Returns:
            list(fluoline)
        """
        if codes is None:
            return []
        elif not hasattr(codes,"__iter__"):
            codes = [codes]

        # line belongs to shell
        # line is selected (directly or under composite)
        # line is not a composite
        use = lambda line,selname,code: line.startswith(shellname) and\
                                        line.startswith(selname) and\
                                        (code<0 or line!=selname)

        lines = [fluoline(cls._all[name]) for code in codes for name in cls._all if use(name,cls.getname(code),code)] # decompose composites

        return list(set(lines)) # remove duplicates


class shell(Hashable):
    _all = {xraylib.__dict__[s]:s.split('_')[0] for s in xraylib.__dict__.keys() if s.endswith("_SHELL")}

    def __init__(self,code,fluolines=None):
        """
        Args:
            code(int): xraylib shell code
            fluolines(Optional(array(int))): xraylib line codes (all lines when None)
        """

        self.code = code
        self.name = self.getname()
        self.setfluo(fluolines=fluolines)

    def _cmpkey(self):
        """For comparing and sorting
        """
        return self.code

    def _stringrepr(self):
        """Unique representation of an instance
        """
        return self.name

    def getname(self):
        return shell._all[self.code]

    def setfluo(self,fluolines=None):
        """
        Args:
            fluolines(Optional(array(int))): xraylib line codes
        """
        self.fluolines = fluoline.factory(fluolines,self.name)

    def fluoyield(self,Z):
        return xraylib.FluorYield(Z,self.code)

    def radrate(self,Z):
        if len(self.fluolines)==0:
            return [1] # ALL lines
        else:
            return [l.radrate(Z) for l in self.fluolines]

    def fluofrac(self,Z):
        return self.fluoyield(Z)*sum(self.radrate(Z))

    def edgenergy(self,Z):
        return xraylib.EdgeEnergy(Z,self.code)


class element(Hashable):
    """Interface to chemical elements
    """

    def __init__(self,symb):
        """symb is a string or an integer
        """
        
        # Atomic number
        if isinstance(symb,str):
            if symb.isdigit():
                self.Z = int(symb)
                self.name = xraylib.AtomicNumberToSymbol(self.Z)
            else:
                self.Z = xraylib.SymbolToAtomicNumber(symb.title())
                self.name = symb
        elif isinstance(symb,int):
            self.Z = symb
            self.name = xraylib.AtomicNumberToSymbol(self.Z)
        else:
            raise ValueError("Element symbol Z.")
        
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

    def markasabsorber(self,symb,shells=None,fluolines=None):
        """Marking an element's shells and lines has following effect:
            - replace cross-section data with simulations for the selected shells
            - partial absorption cross-section is not zero
            - when no fluolines are given: all lines for each shell are taken into account

        Args:
            symb(str): element symbol
            shells(Optional(array(int))): xraylib shell codes
            fluolines(Optional(array(int))): xraylib line codes
        """
        if self.name==symb:
            # Shell names
            if shells is None:
                self.shells = []
            elif not hasattr(shells,"__iter__"):
                self.shells = [shell(shells,fluolines=fluolines)]
            else:
                self.shells = list(set([shell(s,fluolines=fluolines) for s in shells]))

    def unmarkabsorber(self):
        self.shells = []

    def isabsorber(self):
        return len(self.shells)!=0

    def get_pure_density(self):
        return xraylib.ElementDensity(self.Z)

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
        bnum = not hasattr(E,"__iter__")
        if bnum:
            E = [E]
        cs = np.empty(len(E),dtype=np.float64)

        # Total
        for i in range(len(E)):
            cs[i] = xraylib.CS_Total_Kissel(self.Z,E[i])

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
        bnum = not hasattr(E,"__iter__")
        if bnum:
            E = [E]
        cs = np.empty(len(E),dtype=np.float64)

        # Total
        for i in range(len(E)):
            cs[i] = xraylib.CS_Photo_Total(self.Z,E[i])

        # Replace part by simulation
        cs = self._replace_partial_mass_abs_coeff(cs,E,environ=environ,decimals=decimals,refresh=refresh)

        if bnum:
            return cs[0]
        else:
            return cs

    def _replace_partial_mass_abs_coeff(self,cs,E,environ=None,decimals=6,refresh=False):
        """
        """
        if environ is not None and self.isabsorber():
            # Subtract partial cross-sections (summed over selected shells)
            cs -= sum(self._CS_Photo_Partial_DB(E).values())

            ind = np.argwhere(ret<0)
            if len(ind)>0:
                ind2 = np.argwhere(ret>=0)
                f = interpolate.interp1d(E[ind2].flatten(),ret[ind2].flatten(),bounds_error=False)
                ret[ind] = f(E[ind])

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
        bnum = not hasattr(E,"__iter__")
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

        cs = sum(cs.values())
        if bnum:
            return cs[0]
        else:
            return cs

    def xrf_cross_section(self,E,**kwargs):
        return self._xrf_cross_section(E,decomposed=False,**kwargs)

    def xrf_cross_section_decomposed(self,E,**kwargs):
        return self._xrf_cross_section(E,decomposed=True,**kwargs)

    def _xrf_cross_section(self,E,environ=None,decimals=6,refresh=False,decomposed=False):
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
        bnum = not hasattr(E,"__iter__")
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
            cs = sum([c*s.fluofrac(self.Z) for s,c in cs.items()])
            if bnum:
                return cs[0]
            else:
                return cs

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
        if hasattr(E,"__iter__"):
            ret = np.empty(len(E),dtype=np.float64)

            for i in range(len(E)):
                ret[i] = xraylib.CS_Rayl(self.Z,E[i])+xraylib.CS_Compt(self.Z,E[i])
        else:
            ret = xraylib.CS_Rayl(self.Z,E)+xraylib.CS_Compt(self.Z,E)

        return ret

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
        if hasattr(E,"__iter__"):
            ret = np.empty(len(E),dtype=np.float64)

            for i in range(len(E)):
                ret[i] = xraylib.CS_Rayl(self.Z,E[i])
        else:
            ret = xraylib.CS_Rayl(self.Z,E)

        return ret

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
        if hasattr(E,"__iter__"):
            ret = np.empty(len(E),dtype=np.float64)

            for i in range(len(E)):
                ret[i] = xraylib.CS_Compt(self.Z,E[i])
        else:
            ret = xraylib.CS_Compt(self.Z,E)

        return ret

    def _get_cif_name(self,name):
        """Get file from the database if it doesn't exist
        """
        if os.path.isfile(name):
            return name
        else:
            f = os.path.join(os.path.dirname(os.path.abspath(__file__)),"cif",os.path.splitext(os.path.basename(name))[0]+".cif")
            if os.path.isfile(f):
                return f
            else:
                return None

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
            dict: shell:tau(E,shell)
        """
        cs = {}
        for s in self.shells:
            cs[s] = np.empty(len(E),dtype=np.float64)
            for i in range(len(E)):
                cs[s][i] = xraylib.CS_Photo_Partial(self.Z,s.code,E[i])

        return cs

    def _CS_Photo_Partial_SIM(self,E,environ,decimals=6,refresh=False,fluo=True):
        """Calculate the partial photoionization cross section with fdmnes (E in keV).

        Args:
            E(array): energy (keV)
            environ(dict): chemical environment of this element
            decimals(Optional(num)): precision of energy in keV
            refresh(Optional(bool)): force re-simulation if used

        Returns
            dict: shell:tau(E,shell)
        """

        # Initialize simulation
        sim = fdmnes.fdmnes(environ.ciffile, resonant=self.name)

        # Energy range
        sim.P.Energpho = False # relative energies as output

        # Simultation settings
        #sim.P.Radius = 3.5 # Radius of the cluster for calculation
        sim.P.Radius = 5.
        sim.P.Rpotmax = sim.P.Radius + 5 # Radius of the cluster for potential calculation
        sim.P.Quadrupole = True # multipole approximation
        sim.P.Green = False # MS instead of FDM (faster but less accurate)
        sim.P.TDDFT = True # multi electron correction
        sim.P.Convolution = True # 
        sim.P.Density = False # save density of states

        # Do simulation
        cs = {}

        for s in self.shells:
            cs[s] = np.empty(len(E),dtype=np.float64)

            # Select edge
            sim.P.Edge = s.name
            filebase = os.path.splitext(os.path.basename(environ.ciffile))[0]+"_"+self.name+"_"+sim.P.Edge

            # Relative energy range in eV
            Range = self._get_fdmnes_energyrange(E,s.edgeenergy(self.Z),decimals=decimals)
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
            data[:,0] += s.edgeenergy(self.Z) # rel -> abs
            data[:,1] *= xraylib.AVOGNUM*1E6/self.MM # Absorption cross section (Mbarn) -> mass absorption coefficient (cm^2/g)

            # Element multiplicity in the unit cell
            if "nmult" in config:
                nmult = config["nmult"]
            else:
                # TODO: cctbx is too heavy (information is also in the bav file)
                nmult = self._get_multiplicity(environ.structure)
                #nmult = np.round(data[-1,1]/xraylib.CS_Photo_Partial(self.Z,shell,E[-1])) # not precise enough
                config["nmult"] = nmult
            data[:,1] /= nmult

            # Energy interpolation and keep
            f = interpolate.interp1d(data[:,0],data[:,1],bounds_error=False)
            cs[s] = f(E)
        
            # Write configuration file
            with open(fcfg,'w') as f:
                json.dump(config,f,indent=2)

        return ret

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

