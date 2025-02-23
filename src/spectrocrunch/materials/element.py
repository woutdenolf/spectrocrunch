import os
import tempfile
import time
import json
import re

from . import elementbase
from ..patch import xraylib
from ..patch.pint import ureg
from ..utils import instance
from ..math import lazy
from . import xrayspectrum
from . import csutils

import numpy as np
from scipy import interpolate
import warnings
import fisx

try:
    import fdmnes
except ImportError:
    fdmnes = None
    warnings.warn("fdmnes is not installed", ImportWarning)
try:
    import iotbx.cif as iotbxcif
except ImportError:
    iotbxcif = None
    warnings.warn("cctbx is not installed", ImportWarning)


def elementParse(symb):
    if instance.isstring(symb):
        if symb.isdigit():
            Z = int(symb)
            name = xraylib.AtomicNumberToSymbol(Z)
        else:
            Z = xraylib.SymbolToAtomicNumber(str(symb))
            name = symb
    elif instance.isinteger(symb):
        Z = symb
        name = xraylib.AtomicNumberToSymbol(Z)
    elif isinstance(symb, Element):
        Z, name = symb.Z, symb.name
    else:
        raise ValueError("Unknown element symbol or number")
    return Z, name


def elementZ(symb):
    if instance.isstring(symb):
        if symb.isdigit():
            Z = int(symb)
        else:
            Z = xraylib.SymbolToAtomicNumber(symb)
    elif instance.isinteger:
        Z = symb
    elif isinstance(symb, Element):
        Z = symb.Z
    else:
        raise ValueError("Unknown element symbol or number")
    return Z


def elementSymbol(symb):
    if instance.string(symb):
        if symb.isdigit():
            name = xraylib.AtomicNumberToSymbol(int(symb))
        else:
            name = symb
    elif instance.isinteger:
        name = xraylib.AtomicNumberToSymbol(symb)
    elif isinstance(symb, Element):
        name = symb.name
    else:
        raise ValueError("Unknown element symbol or number")
    return name


class Element(elementbase.ElementBase):
    """Interface to chemical elements"""

    def __init__(self, symb):
        """
        Args:
            symb(int or str or Element)
        """
        self.Z, self.name = elementParse(symb)
        self.shells = []
        self.isscatterer = True

    @property
    def Z(self):
        return self._Z

    @Z.setter
    def Z(self, value):
        self._Z, self._name = elementParse(value)
        self.MM = xraylib.AtomicWeight(self._Z)

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        self.Z = value

    def __getstate__(self):
        return {"Z": self.Z, "shells": self.shells, "isscatterer": self.isscatterer}

    def __setstate__(self, state):
        self.Z = state["Z"]
        self.shells = state["shells"]
        self.isscatterer = state["isscatterer"]

    def _sortkey(self, other):
        return -self.Z

    def _cmpkey(self, other):
        if instance.isnumber(other):
            return self.Z
        else:
            return repr(self)

    @property
    def _repr(self):
        return self.name

    def shellfactory(self, emin=None, emax=None):
        return xrayspectrum.Shell.factory(energybounds=[self.Z, emin, emax])

    def pymcashellfactory(self, emin=None, emax=None):
        return xrayspectrum.Shell.pymcafactory(energybounds=[self.Z, emin, emax])

    def markabsorber(self, symb=None, shells=None, fluolines=None, energybounds=None):
        """
        This is needs for the partial absorption and XRF cross-section (zero otherwise).

        Args:
            symb(str): element symbol
            shells(Optional(array(int))): restrict shells
            fluolines(Optional(array(int))): restrict lines (`None`: all, `[]`:explicite all)
            energybounds(Optional(2-tuple)): restrict shells based on energy (only when `shells` is `None`)
        """
        if symb is None:
            mark = True
        else:
            Z, name = elementParse(symb)
            mark = self.Z == Z
        if mark:
            if shells is None:
                if energybounds is not None:
                    shells = self.shellfactory(
                        emin=energybounds[0], emax=energybounds[1]
                    )
            if shells is None:
                self.shells = xrayspectrum.Shell.all_shells(fluolines=fluolines)
            else:
                shells = xrayspectrum.Shell.expandall(shells)

                def func(shell):
                    if isinstance(shell, xrayspectrum.Shell):
                        return shell
                    else:
                        return xrayspectrum.Shell(shell, fluolines=fluolines)

                if instance.isiterable(shells):
                    self.shells = [func(shell) for shell in shells]
                else:
                    self.shells = [func(shells)]

    def edge_energies(self):
        return [shell.edgeenergy(self.Z) for shell in self.shells]

    def unmarkabsorber(self):
        self.shells = []

    @property
    def isabsorber(self):
        return bool(self.shells)

    def markscatterer(self, name=None):
        if name is None:
            self.isscatterer = True
        else:
            self.isscatterer = self == name

    def unmarkscatterer(self):
        self.isscatterer = False

    def markinfo(self):
        if self.isabsorber:
            yield "{} ionized shells:".format(self.name)
            for shell in self.shells:
                for s in shell.markinfo():
                    yield " {}".format(s)
        else:
            yield "{} (no ionized shells)".format(self.name)

    @property
    def density(self):
        return xraylib.ElementDensity(self.Z)

    def molarmass(self):
        return self.MM

    def molarmasseff(self):
        return self.MM

    def massfractions(self):
        return {self: 1.0}

    def molefractions(self):
        return self.equivalents()

    def equivalents(self):
        return {self: 1.0}

    def elemental_molefractions(self):
        return self.molefractions()

    def elemental_equivalents(self):
        return self.equivalents()

    def elemental_massfractions(self):
        return self.massfractions()

    @property
    def nelements(self):
        return 1

    @property
    def elements(self):
        return [self]

    @property
    def parts(self):
        return {self: 1}

    def _xraylib_method(self, method, E):
        return csutils.eval(method, self.Z, E, applypost=False)

    def _xraylib_method_full(self, method, E):
        return csutils.eval(method, self.Z, E, applypost=True)

    def mass_att_coeff(self, E, environ=None, decimals=6, refresh=False, **kwargs):
        """
        Mass attenuation coefficient (cm^2/g, E in keV). Use for transmission XAS.

        Args:
            E(num or array-like): energy (keV)
            environ(dict): chemical environment of this element
            decimals(Optional(num)): precision of energy in keV
            refresh(Optional(bool)): force re-simulation if used

        Returns:
            num or np.array
        """
        # Total
        cs, func = self._xraylib_method("CS_Total_Kissel", E)
        # Replace part by simulation
        cs = self._replace_partial_mass_abs_coeff(
            cs, E, environ=environ, decimals=decimals, refresh=refresh
        )
        return func(cs)

    def mass_abs_coeff(self, E, environ=None, decimals=6, refresh=False, **kwargs):
        """
        Mass absorption coefficient (cm^2/g, E in keV).

        Args:
            E(num or array-like): energy (keV)
            environ(dict): chemical environment of this element
            decimals(Optional(num)): precision of energy in keV
            refresh(Optional(bool)): force re-simulation if used

        Returns:
            num or np.array
        """
        # Total
        cs, func = self._xraylib_method("CS_Photo_Total", E)
        # Replace part by simulation
        cs = self._replace_partial_mass_abs_coeff(
            cs, E, environ=environ, decimals=decimals, refresh=refresh
        )
        return func(cs)

    def _replace_partial_mass_abs_coeff(
        self, cs, E, environ=None, decimals=6, refresh=False, **kwargs
    ):
        """"""
        if environ is not None:
            E = np.atleast_1d(E)
            # Subtract partial cross-sections (summed over selected shells)
            cs -= sum(self._CS_Photo_Partial_DB(E).values())
            ind = np.argwhere(cs < 0)
            if len(ind) > 0:
                ind2 = np.argwhere(cs >= 0)
                f = interpolate.interp1d(
                    E[ind2].flatten(), cs[ind2].flatten(), bounds_error=False
                )
                cs[ind] = f(E[ind])
            # Add partial cross-sections (summed over selected shells)
            cs += sum(
                self._CS_Photo_Partial_SIM(
                    E, environ, decimals=decimals, refresh=refresh
                ).values()
            )
        return cs

    def partial_mass_abs_coeff(
        self, E, environ=None, decimals=6, refresh=False, **kwargs
    ):
        """
        Mass absorption coefficient for the selected shells (cm^2/g, E in keV).

        Args:
            E(num or array-like): energy (keV)
            environ(dict): chemical environment of this element
            decimals(Optional(num)): precision of energy in keV
            refresh(Optional(bool)): force re-simulation if used

        Returns:
            num or array: sum_S[tau(E,S)]
        """
        E, func = instance.asarrayf(E)

        if not self.isabsorber:
            return func(np.zeros(len(E), dtype=np.float64))

        if environ is None:
            cs = self._CS_Photo_Partial_DB(E)
        else:
            cs = self._CS_Photo_Partial_SIM(
                E, environ, decimals=decimals, refresh=refresh
            )

        # sum over selected shells
        cs = sum(cs.values())
        return func(cs)

    def fluorescence_cross_section(self, E, environ=None, decomposed=False, **kwargs):
        r"""
        XRF cross section for the selected shells and lines (cm^2/g, E in keV). Use for fluorescence XAS.

        Args:
            E(num or array-like): energy (keV)
            environ(dict): chemical environment of this element
            decomposed(Optional(bool)): output as dataframe
            \**kwargs: see _CS_Photo_Partial_SIM

        Returns:
            num or np.array: sum_{S}[tau(E,S)*fluoyield(S)*sum_{L}[radrate(S,L)]]
            dict: S:tau(E,S)
        """
        E, func = instance.asarrayf(E)
        if self.isabsorber:
            if environ is None:
                cs = self._CS_Photo_Partial_DB(E)
            else:
                cs = self._CS_Photo_Partial_SIM(E, environ, **kwargs)
            if not decomposed:
                # weighted sum over selected shells, weight is the total fluoyield of the selected lines belonging to the shell
                cs = sum(
                    [
                        shellcs * shell.partial_fluoyield(self.Z)
                        for shell, shellcs in cs.items()
                    ]
                )
                cs = func(cs)
        else:
            if decomposed:
                cs = {}
            else:
                cs = func(np.zeros(len(E), dtype=np.float64))
        return cs

    def fluorescence_cross_section_lines(self, E, **kwargs):
        r"""
        XRF cross sections per line (cm^2/g, E in keV). Use for XRF.

        Args:
            E(num or array-like): energy (keV)
            \**kwargs: see fluorescence_cross_section

        Returns:
            dict: L:tau(E,S)*fluoyield(S)*radrate(S,L)
        """
        if self.isabsorber:
            # Get the shell ionization cross section
            kwargs["decomposed"] = True
            cs = self.fluorescence_cross_section(E, **kwargs)
            cs = {
                xrayspectrum.FluoZLine(self, line): shellcs * fluoyield
                for shell, shellcs in cs.items()
                for line, fluoyield in shell.partial_fluoyield(
                    self.Z, decomposed=True
                ).items()
            }
        else:
            cs = {}
        return cs

    def diff_fluorescence_cross_section(self, E, **kwargs):
        """Differential XRF cross sections per line (cm^2/g/srad, E in keV). Use for XRF."""
        cs = self.fluorescence_cross_section_lines(E, **kwargs)
        m = 4 * np.pi
        return {k: v / m for k, v in cs.items()}

    def scattering_cross_section(self, E, **kwargs):
        """Scattering cross section (cm^2/g, E in keV).

        Args:
            E(num or array-like): energy (keV)

        Returns:
            num or np.array
        """
        return self._xraylib_method_full("CS_Rayl", E) + self._xraylib_method_full(
            "CS_Compt", E
        )

    def rayleigh_cross_section(self, E, **kwargs):
        """Rayleigh cross section (cm^2/g, E in keV).

        Args:
            E(num or array-like): energy (keV)

        Returns:
            num or np.array
        """
        return self._xraylib_method_full("CS_Rayl", E)

    def compton_cross_section(self, E, **kwargs):
        """Compton cross section (cm^2/g, E in keV).

        Args:
            E(num or array-like): energy (keV)

        Returns:
            num or np.array
        """
        return self._xraylib_method_full("CS_Compt", E)

    def scatfact_classic_real(self, E, theta=None, **kwargs):
        """Real part of atomic form factor

        Args:
            E(num or array-like): energy (keV)
            theta(Optional(num or array-like)): scattering angle (rad)

        Returns:
            num or np.array
        """
        if theta is None:
            return self.Z
        else:
            # q = Q/4.pi
            q = (
                np.sin(theta / 2)
                / ureg.Quantity(E, "keV").to("angstrom", "spectroscopy").magnitude
            )
            return self._xraylib_method_full("FF_Rayl", q)

    def scatfact_real(self, E, theta=None, **kwargs):
        """Real part of atomic form factor

        Args:
            E(num or array-like): energy (keV)
            theta(Optional(num or array-like)): scattering angle (degrees)

        Returns:
            num or np.array
        """
        return self.scatfact_classic_real(E, theta=theta) + self._xraylib_method_full(
            "Fi", E
        )

    def scatfact_imag(self, E, **kwargs):
        """Imaginary part of atomic form factor

        Args:
            E(num or array-like): energy (keV)

        Returns:
            num or np.array
        """
        return self._xraylib_method_full("Fii", E)

    def diff_rayleigh_cross_section(self, E, source=None, **kwargs):
        """Differential Rayleigh cross section (cm^2/g/srad, E in keV).

        Args:
            E(num or array-like): energy (keV)
            source(Source): X-ray source

        Returns:
            callable: (azimuth,polar)
        """
        # muR(cm²/g) = NA(atom/mol)/MM(g/mol) . x(cm²/atom)
        # x(cm²/atom) = int_phi[int_theta [ f²(e/atom).thomson(cm²/e/srad) . sin(theta) dtheta] dphi]
        # thomson(cm²/e/srad) = r_e²(cm²/e).K_thomson(azimuth,polar)(1/srad)
        # mudiffR(cm²/g/srad) = r_e²(cm²/e).NA(atom/mol)/MM(g/mol).K_thomson(azimuth,polar)(1/srad).f²(e/atom)

        c = (
            (
                ureg.classical_electron_radius**2
                * ureg.particles_per_mol
                / ureg.Quantity(self.MM, "g/mol")
            )
            .to("cm^2/g")
            .magnitude
        )
        wl = ureg.Quantity(E, "keV").to("angstrom", "spectroscopy").magnitude
        K = source.thomson_K

        def func(azimuth, polar):
            return (
                c
                * K(azimuth, polar)
                * self._xraylib_method_full("FF_Rayl", np.sin(polar / 2.0) / wl) ** 2
            )

        return lazy.Function(func, name="diff_el(phi,theta)")

    def diff_compton_cross_section(self, E, source=None, **kwargs):
        """Differential Compton cross section (cm^2/g/srad, E in keV).

        Args:
            E(num or array-like): energy (keV)
            source(Source): X-ray source
            decomposed(Optional(bool)): not used

        Returns:
            callable: (theta,phi)
        """
        # muC(cm²/g) = NA(atom/mol)/MM(g/mol) . x(cm²/atom)
        # x(cm²/atom) = int_phi[int_theta [ S(e/atom).KN(cm²/e/srad) . sin(theta) dtheta] dphi]
        # KN(cm²/e/srad) = r_e²(cm²/e).K_compton(azimuth,polar)(1/srad)
        # mudiffC(cm²/g/srad) = r_e²(cm²/e).NA(atom/mol)/MM(g/mol).K_compton(azimuth,polar)(1/srad).S(e/atom)
        c = (
            (
                ureg.classical_electron_radius**2
                * ureg.particles_per_mol
                / ureg.Quantity(self.MM, "g/mol")
            )
            .to("cm^2/g")
            .magnitude
        )
        wl = ureg.Quantity(E, "keV").to("angstrom", "spectroscopy").magnitude
        K = source.compton_K(E)

        def func(azimuth, polar):
            return (
                c
                * K(azimuth, polar)
                * self._xraylib_method_full("SF_Compt", np.sin(polar / 2.0) / wl)
            )

        return lazy.Function(func, name="diff_inel(phi,theta)")

    def _get_multiplicity(self, struct):
        scat = struct.scatterers()
        ret = 0.0
        for s in scat:
            if s.scattering_type == self.name:
                ret += s.occupancy * s.multiplicity()
        return ret

    def _CS_Photo_Partial_DB(self, E):
        """
        Get the partial photoionization cross section for the selected shells from xraylib (E in keV).

        Args:
            E(array): energy (keV)

        Returns
            dict: xrayspectrum.Shell:tau(E,xrayspectrum.Shell)
        """
        cs = {}
        for shell in self.shells:
            cs[shell] = np.vectorize(lambda en: shell.partial_photo(self.Z, en))(E)
        return cs

    def _CS_Photo_Partial_SIM(self, E, environ, decimals=6, refresh=False, fluo=True):
        """
        Calculate the partial photoionization cross section for the selected shells from fdmnes (E in keV).

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
        sim.P.Energpho = False  # relative energies as output

        # Simultation settings
        # sim.P.Radius = 3.5 # Radius of the cluster for calculation
        sim.P.Radius = 7.0
        sim.P.Rpotmax = (
            sim.P.Radius + 5
        )  # Radius of the cluster for potential calculation
        sim.P.Quadrupole = True  # multipole approximation
        sim.P.Green = False  # MS instead of FDM (faster but less accurate)
        sim.P.TDDFT = True  # multi electron correction
        sim.P.Convolution = True
        sim.P.Density = False  # save density of states

        # Do simulation
        cs = {}

        for shell in self.shells:
            cs[shell] = np.empty(len(E), dtype=np.float64)

            # Select edge
            sim.P.Edge = shell.name
            filebase = "_".join(
                (
                    os.path.splitext(os.path.basename(environ.ciffile))[0],
                    self.name,
                    sim.P.Edge,
                )
            )

            # Relative energy range in eV
            Range = self._get_fdmnes_energyrange(
                E, shell.edgeenergy(self.Z), decimals=decimals
            )
            sim.P.Range = tuple(Range)

            # Read previous configuration file
            fcfg = os.path.join(tempfile.gettempdir(), filebase + ".json")
            if os.path.isfile(fcfg) and not refresh:
                with open(fcfg, "r") as f:
                    config = json.load(f)
            else:
                config = {}

            # Write simulation file
            finput = os.path.join(tempfile.gettempdir(), filebase + ".txt")
            sim.WriteInputFile(finput, overwrite=True)

            # Perform simulation if not already done
            if "Range" in config:
                Range_prev = np.array(config["Range"])
            else:
                Range_prev = Range * 0
            config["Range"] = sim.P.Range

            if len(Range_prev) == len(Range):
                brun = not np.allclose(Range_prev, Range)
            else:
                brun = True
            if brun:
                sim.Run(wait=False)
                while True:
                    NumRunning = sum([not sim.Status(i) for i in range(len(sim.proc))])
                    if NumRunning == 0:
                        break
                    time.sleep(5)

            # Get data
            data = sim.get_XANES(conv=True)

            # Convert data to absorbance as a function of energy in keV
            data[:, 0] /= 1000.0  # ev -> keV
            data[:, 0] += shell.edgeenergy(self.Z)  # rel -> abs
            data[:, 1] *= xraylib.AVOGNUM * 1e6 / self.MM  # Mbarn/atom -> cm^2/g

            # Element multiplicity in the unit cell
            if "nmult" in config:
                nmult = config["nmult"]
            else:
                # TODO: cctbx is too heavy (information is also in the bav file)
                nmult = self._get_multiplicity(environ.structure)
                # nmult = np.round(data[-1,1]/xraylib.CS_Photo_Partial(self.Z,Shell,E[-1])) # not precise enough
                config["nmult"] = nmult
            data[:, 1] /= nmult

            # Energy interpolation and keep
            f = interpolate.interp1d(data[:, 0], data[:, 1], bounds_error=False)
            cs[shell] = f(E)

            # Write configuration file
            with open(fcfg, "w") as f:
                json.dump(config, f, indent=2)

        return cs

    def get_energy(self, energyrange, defaultinc=1):
        """Get absolute energies (keV) from a relative energy range (eV)"""
        nshells = len(self.shells)
        if nshells == 0:
            return None
        elif nshells == 1:
            return (
                self.shells[0].edgeenergy(self.Z)
                + self._get_fdmnes_energies(energyrange) / 1000.0
            )
        else:
            return self._get_fdmnes_energies(
                self._get_absolute_energyrange(energyrange, defaultinc=defaultinc)
            )

    def _get_fdmnes_energies(self, energyrange):
        """Calculate energies based on boundaries and step sizes:
        energyrange = [E0,step0,E1,step1,E2]
        """

        # Number of steps in each region
        nblocks = len(energyrange) // 2
        nsteps = np.empty(nblocks, dtype=int)
        e = energyrange[0]
        for i in range(nblocks):
            b = e
            e = energyrange[2 * i + 2]
            inc = energyrange[2 * i + 1]
            nsteps[i] = np.ceil((e - b) / inc)
            e = b + nsteps[i] * inc
            if i == nblocks - 1:
                if e > energyrange[-1]:
                    nsteps[i] -= 1

        ret = np.empty(nsteps.sum() + 1, dtype=energyrange.dtype)
        ret[0] = energyrange[0]
        off = 1
        for i in range(nblocks):
            inc = energyrange[2 * i + 1]
            ret[off : off + nsteps[i]] = (
                ret[off - 1] + np.arange(1, nsteps[i] + 1) * inc
            )
            off += nsteps[i]

        return ret

    def _get_fdmnes_energyrange(self, Eabs, edgeenergy, decimals=6):
        """Calculate energies based on boundaries and step sizes:
         energyrange = [E0,step0,E1,step1,E2] (eV, relative to the edge)
        Eabs in keV
        """

        E = (Eabs - edgeenergy) * 1000  # absolute (keV) to relative (eV)

        dE = np.around(E[1:] - E[0:-1], decimals=decimals - 3)
        ind = np.argwhere(dE[1:] - dE[0:-1])
        nblocks = len(ind) + 1

        energyrange = np.empty(2 * nblocks + 1, dtype=E.dtype)
        energyrange[0] = E[0]
        energyrange[-1] = E[-1]

        if nblocks == 1:
            energyrange[1] = dE[0]
        else:
            inddest = np.arange(1, 2 * nblocks - 2, 2)
            energyrange[inddest] = dE[ind]
            energyrange[inddest + 1] = E[ind + 1]
            energyrange[-2] = dE[ind[-1] + 1]

        # enlarge the fdmnes range a bit
        add = 5  # eV
        energyrange[0] -= np.ceil(add / energyrange[1]) * energyrange[1]
        energyrange[-1] += np.ceil(add / energyrange[-2]) * energyrange[-2]

        return energyrange

    def _get_absolute_energyrange(self, energyrange, defaultinc=1):
        """Convert relative energy range (eV) to absolute energy range (keV)"""

        # Boundaries
        indE = np.arange(0, len(energyrange), 2)
        nshells = len(self.shells)
        nbounds = len(indE)
        energies = np.empty(nshells * nbounds, dtype=energyrange.dtype)
        for i in range(nshells):
            energies[i * nbounds : (i + 1) * nbounds] = (
                self.shells[i].edgeenergy(self.Z) + energyrange[indE] / 1000.0
            )

        # Put unique boundaries in new range array
        E = np.unique(energies)
        newnblocks = len(E) - 1
        newrange = np.empty(2 * newnblocks + 1, dtype=energyrange.dtype)
        newindE = np.arange(0, len(newrange), 2)
        newindD = np.arange(1, 2 * newnblocks, 2)
        newrange[newindE] = E
        newrange[newindD] = defaultinc / 1000.0

        # Determine step sizes
        for j in range(newnblocks):
            m = newrange[2 * j] + (newrange[2 * j + 2] - newrange[2 * j]) / 2
            inc = 0

            for i in range(nshells):
                b = i * nbounds
                e = b + nbounds - 1
                if m > energies[b] and m < energies[e]:
                    ubound = 2 * np.argmax((energies[b : e + 1] - m) > 0)
                    d = energyrange[ubound - 1] / 1000.0
                    if inc == 0:
                        inc = d
                    else:
                        inc = min(inc, d)

            if inc != 0:
                newrange[2 * j + 1] = inc

        return newrange

    def topymca(self, cfg, defaultthickness=1e-4):
        r = self.massfractions()
        massfractions = list(r.values())
        names = ["{}1".format(e) for e in r]
        matname = self.pymcaname
        cfg["materials"][matname] = {
            "Comment": self.pymcacomment,
            "CompoundFraction": massfractions,
            "Thickness": defaultthickness,
            "Density": self.density,
            "CompoundList": names,
        }
        return matname

    def tofisx(self, cfg, defaultthickness=1e-4):
        r = self.massfractions()
        massfractions = list(r.values())
        names = ["{}1".format(e) for e in r]
        matname = self.pymcaname
        o = fisx.Material(matname, self.density, defaultthickness, self.pymcacomment)
        o.setCompositionFromLists(names, massfractions)
        cfg.addMaterial(o, errorOnReplace=False)
        return matname

    @classmethod
    def fluozgroup(cls, symb):
        ele, group = re.split("[-_ ]+", symb)
        return xrayspectrum.FluoZGroup(cls(ele), group)

    def fluolines(self, lines):
        if instance.isstring(lines):
            return xrayspectrum.FluoZLine(self, xrayspectrum.FluoLine(lines))
        else:
            return [
                xrayspectrum.FluoZLine(self, xrayspectrum.FluoLine(line))
                for line in lines
            ]

    def tocompound(self, name=None):
        return self
