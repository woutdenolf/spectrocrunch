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

try:
    import iotbx.cif
except:
    pass

import fdmnes
import fisx

from .. import xraylib
from .. import ureg
from ..common import hashable
from ..common import instance
from . import xrayspectrum


def elementParse(symb):
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

def elementZ(symb):
    if isinstance(symb,str):
        if symb.isdigit():
            Z = int(symb)
        else:
            Z = xraylib.SymbolToAtomicNumber(symb.title())
    elif isinstance(symb,numbers.Integral):
        Z = symb
    elif isinstance(symb,Element):
        Z = symb.Z
    else:
        raise ValueError("Unknown element symbol or number")
    return Z

def elementSymbol(symb):
    if isinstance(symb,str):
        if symb.isdigit():
            name = xraylib.AtomicNumberToSymbol(int(symb))
        else:
            name = symb
    elif isinstance(symb,numbers.Integral):
        name = xraylib.AtomicNumberToSymbol(symb)
    elif isinstance(symb,Element):
        name = symb.name
    else:
        raise ValueError("Unknown element symbol or number")
    return name  

class Element(hashable.Hashable):
    """Interface to chemical elements
    """

    def __init__(self,symb):
        """
        Args:
            symb(int or str or Element)
        """
        
        # Atomic number
        self.Z,self.name = elementParse(symb)
        
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

    def shellfactory(self,emin=None,emax=None):
        return xrayspectrum.Shell.factory(energybounds=[self.Z,emin,emax])

    def pymcashellfactory(self,emin=None,emax=None):
        return xrayspectrum.Shell.pymcafactory(energybounds=[self.Z,emin,emax])

    def markabsorber(self,symb=None,shells=None,fluolines=None,energybounds=None):
        """Marking an element's shells and lines has following effect:
            - partial absorption cross-section is not zero
            - when no fluolines are given: all lines for each shell are taken into account

        Args:
            symb(str): element symbol
            shells(Optional(array(int))): None -> all
            fluolines(Optional(array(int))): None -> all, [] -> explicite all
        """

        if symb is None:
            mark = True
        else:
            Z,name = elementParse(symb)
            mark = self.Z==Z

        if mark:
            if shells is None:
                if energybounds is not None:
                    shells = self.shellfactory(emin=energybounds[0],emax=energybounds[1])
        
            if shells is None:
                self.shells = xrayspectrum.Shell.all_shells(fluolines=fluolines)
            else:
                f = lambda shell: shell if isinstance(shell,xrayspectrum.Shell) else xrayspectrum.Shell(shell,fluolines=fluolines)
                if instance.isiterable(shells):
                    self.shells = [f(shell) for shell in shells]
                else:
                    self.shells = [f(shells)]

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

    @property
    def ncompounds(self):
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

    def mass_att_coeff(self,E,environ=None,decimals=6,refresh=False,**kwargs):
        """Mass attenuation coefficient (cm^2/g, E in keV). Use for transmission XAS.

        Args:
            E(num or array-like): energy (keV)
            environ(dict): chemical environment of this element
            decimals(Optional(num)): precision of energy in keV
            refresh(Optional(bool)): force re-simulation if used

        Returns:
            num or np.array
        """
        
        E,func = instance.asarrayf(E)
        cs = np.empty(len(E),dtype=np.float64)

        # Total
        for i,en in enumerate(E):
            cs[i] = xraylib.CS_Total_Kissel(self.Z,np.float64(en))

        # Replace part by simulation
        cs = self._replace_partial_mass_abs_coeff(cs,E,environ=environ,decimals=decimals,refresh=refresh)

        return func(cs)

    def mass_abs_coeff(self,E,environ=None,decimals=6,refresh=False,**kwargs):
        """Mass absorption coefficient (cm^2/g, E in keV).

        Args:
            E(num or array-like): energy (keV)
            environ(dict): chemical environment of this element
            decimals(Optional(num)): precision of energy in keV
            refresh(Optional(bool)): force re-simulation if used

        Returns:
            num or np.array
        """
        E,func = instance.asarrayf(E)

        # Total
        cs = np.vectorize(lambda en:xraylib.CS_Photo_Total(self.Z,np.float64(en)))(E)

        # Replace part by simulation
        cs = self._replace_partial_mass_abs_coeff(cs,E,environ=environ,decimals=decimals,refresh=refresh)

        return func(cs)

    def _replace_partial_mass_abs_coeff(self,cs,E,environ=None,decimals=6,refresh=False,**kwargs):
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

    def partial_mass_abs_coeff(self,E,environ=None,decimals=6,refresh=False,**kwargs):
        """Mass absorption coefficient for the selected shells (cm^2/g, E in keV).

        Args:
            E(num or array-like): energy (keV)
            environ(dict): chemical environment of this element
            decimals(Optional(num)): precision of energy in keV
            refresh(Optional(bool)): force re-simulation if used

        Returns:
            num or array: sum_S[tau(E,S)]
        """
        E,func = instance.asarrayf(E)

        if not self.isabsorber():
            return func(np.zeros(len(E),dtype=np.float64))

        if environ is None:
            cs = self._CS_Photo_Partial_DB(E)
        else:
            cs = self._CS_Photo_Partial_SIM(E,environ,decimals=decimals,refresh=refresh)

        # sum over selected shells
        cs = sum(cs.values())
        return func(cs)

    def fluorescence_cross_section(self,E,environ=None,decimals=6,refresh=False,decomposed=False,**kwargs):
        """XRF cross section for the selected shells and lines (cm^2/g, E in keV). Use for fluorescence XAS.

        Args:
            E(num or array-like): energy (keV)
            environ(dict): chemical environment of this element
            decimals(Optional(num)): precision of energy in keV
            refresh(Optional(bool)): force re-simulation if used
            decomposed(Optional(bool)): output as dictionary

        Returns:
            num or np.array: sum_{S}[tau(E,S)*fluoyield(S)*sum_{L}[radrate(S,L)]]
            dict: S:tau(E,S)
        """
        E,func = instance.asarrayf(E)

        if not self.isabsorber():
            return func(np.zeros(len(E),dtype=np.float64))

        if environ is None:
            cs = self._CS_Photo_Partial_DB(E)
        else:
            cs = self._CS_Photo_Partial_SIM(E,environ,decimals=decimals,refresh=refresh)

        if decomposed:
            return {shell:func(shellcs) for shell,shellcs in cs.items()}
        else:
            # sum over selected shells, weighted but the total fluoyield of the selected lines
            cs = sum([shellcs*shell.partial_fluoyield(self.Z) for shell,shellcs in cs.items()])
            return func(cs)

    def fluorescence_cross_section_lines(self,E,environ=None,decimals=6,refresh=False,decomposed=True,**kwargs):
        """XRF cross sections per line (cm^2/g, E in keV). Use for XRF.

        Args:
            E(num or array-like): energy (keV)
            environ(dict): chemical environment of this element
            decimals(Optional(num)): precision of energy in keV
            refresh(Optional(bool)): force re-simulation if used
            decomposed(Optional(bool)): per line or per shell

        Returns:
            dict: {S:tau(E,S)*fluoyield(S)*sum_{L}[radrate(S,L)]}
            dict: {L:tau(E,S)*fluoyield(S)*radrate(S,L)}
        """
        
        # No fluorescence when not an absorber
        if not self.isabsorber():
            return {}

        # Get the shell ionization cross section
        cs = self.fluorescence_cross_section(E,environ=environ,decimals=decimals,refresh=refresh,decomposed=True)
            
        if decomposed:
            # Multiply by fluorescence yield for each line
            cs = {xrayspectrum.FluoZLine(self,line):shellcs*fluoyield\
                       for shell,shellcs in cs.items()\
                       for line,fluoyield in shell.partial_fluoyield(self.Z,decomposed=True).items()}
                       #if line.energy(self.Z)>0
        else:
            cs = {shell:shellcs*shell.partial_fluoyield(self.Z,decomposed=False)\
                        for shell,shellcs in cs.items()}
        return cs

    def xrayspectrum(self,E,weights=None,emin=0,emax=None):
        E = instance.asarray(E)
        if emax is None:
            emax = E[-1]
        self.markabsorber(energybounds=[emin,emax])

        spectrum = xrayspectrum.Spectrum()
        spectrum.density = self.density
        spectrum.update(self.fluorescence_cross_section_lines(E,decomposed=True))
        spectrum[xrayspectrum.RayleighLine(E)] = self.rayleigh_cross_section(E)
        spectrum[xrayspectrum.ComptonLine(E)] = self.compton_cross_section(E)
        spectrum.apply_weights(weights)
        spectrum.xlim = [emin,emax]
        spectrum.title = str(self)
        spectrum.type = spectrum.TYPES.crosssection
        
        return spectrum

    def _xraylib_method(self,method,E):
        method = getattr(xraylib,method)
        E,func = instance.asarrayf(E)
        ret = np.vectorize(lambda en: method(self.Z,np.float64(en)))(E)
        return func(ret)
        
    def scattering_cross_section(self,E,environ=None,decimals=6,refresh=False,**kwargs):
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

    def rayleigh_cross_section(self,E,environ=None,decimals=6,refresh=False,**kwargs):
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

    def compton_cross_section(self,E,environ=None,decimals=6,refresh=False,**kwargs):
        """Compton cross section (cm^2/g, E in keV).

        Args:
            E(num or array-like): energy (keV)
            environ(dict): chemical environment of this element
            decimals(Optional(num)): precision of energy in keV
            refresh(Optional(bool)): force re-simulation if used

        Returns:
            num or np.array
        """
        return self._xraylib_method("CS_Compt",E)

    def scatfact_classic_re(self,E,theta=None,environ=None,decimals=6,refresh=False,**kwargs):
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
            
    def scatfact_re(self,E,theta=None,environ=None,decimals=6,refresh=False,**kwargs):
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
    
    def scatfact_im(self,E,environ=None,decimals=6,refresh=False,**kwargs):
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
            dict: xrayspectrum.Shell:tau(E,xrayspectrum.Shell)
        """
        cs = {}
        for shell in self.shells:
            cs[shell] = np.vectorize(lambda en: shell.partial_photo(self.Z,en))(E)

        return cs

    def _CS_Photo_Partial_SIM(self,E,environ,decimals=6,refresh=False,fluo=True):
        """Calculate the partial photoionization cross section with fdmnes (E in keV).

        Args:
            E(array): energy (keV)
            environ(dict): chemical environment of this element
            decimals(Optional(num)): precision of energy in keV
            refresh(Optional(bool)): force re-simulation if used

        Returns
            dict: xrayspectrum.Shell:tau(E,xrayspectrum.Shell)
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

    @property
    def pymcaname(self):
        return '{}1'.format(self.name)

    def topymca(self,defaultthickness=1e-4):
        value = {'Comment': self.name,
                'CompoundFraction': [1.],
                'Thickness': defaultthickness,
                'Density': self.density,
                'CompoundList': [self.pymcaname]}
        return self.name,value
        
    def tofisx(self):
        o = fisx.Material(self.name, self.density, 1e-10)
        o.setCompositionFromLists([self.pymcaname],[1.])
        return o

    def fisxgroups(self,emin=0,emax=np.inf):
        self.markabsorber(energybounds=[emin,emax])
        el = str(self)
        return {self:self.shells}
        
        
