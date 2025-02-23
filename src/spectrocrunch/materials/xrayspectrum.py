import matplotlib.pyplot as plt
import numpy as np
import numbers
import collections
import collections.abc
import warnings
import scipy.integrate
import copy
from collections import Counter

from ..patch import xraylib
from ..patch.pint import ureg
from ..utils import instance
from ..utils.hashable import CompHashable
from ..utils.copyable import Copyable
from ..utils import listtools
from ..utils.Enum import Enum
from ..utils.roi import mergeroi1d
from ..math import fit1d

# Some xraylib comments:
# CS_FluorLine_Kissel_Cascade(S,X) = FluorYield(S)*RadRate(SX)*PS_full_cascade_kissel
# PS_full_cascade_kissel = CS_Photo_Partial(S) + Cascade(...) + CosKron(...)
# ... -> use PS_full_cascade_kissel,FluorYield,RadRate of lower shells
#
# CS_FluorLine_Kissel_no_Cascade(S,X) = FluorYield(S)*RadRate(SX)*PS_pure_kissel
# PS_full_cascade_kissel = CS_Photo_Partial(S) + CosKron(...)
# ... -> use PS_pure_kissel,FluorYield,RadRate of lower shells


class FluoLine(CompHashable, Copyable):
    @staticmethod
    def getlinename(line):
        # Return IUPAC instead of SIEGBAHN when possible
        if isinstance(line, FluoLine):
            return line.name
        elif isinstance(line, numbers.Number):
            names = xraylib.code_to_line[line]

            n = len(names)
            if n == 1:
                return names[0]
            elif n == 2:
                n0 = len(names[0])
                n1 = len(names[1])
                if n0 == n1:
                    if names[0][1] == "A" or names[0][1] == "B":
                        return names[1]
                    else:
                        return names[0]
                elif n0 > n1:
                    return names[0]
                else:
                    return names[1]
            else:
                return names[np.argmax([len(name) for name in names])]
        else:
            if line not in xraylib.line_to_code.keys():
                raise ValueError("Unknown line name {}".format(line))
            return line

    @staticmethod
    def getlinecode(line):
        if isinstance(line, FluoLine):
            return line.code
        elif isinstance(line, numbers.Number):
            if line not in xraylib.code_to_line.keys():
                raise ValueError("Unknown line code {}".format(line))
            return line
        else:
            return xraylib.line_to_code[line]

    @staticmethod
    def decompose(code):
        """
        Args:
            code(int):
        Return:
            list(int): sub lines
        """
        if code in xraylib.composites:
            return xraylib.composites[code]
        else:
            return [code]

    @classmethod
    def all_lines(cls):
        return list(
            set(
                cls(code)
                for line in range(xraylib.line_max, xraylib.line_min - 1, -1)
                for code in cls.decompose(line)
            )
        )

    @classmethod
    def factory(cls, shells=None, fluolines=None, energybounds=None):
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
                if fluolines == []:
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

            def valid(linecode):
                return any(
                    any(
                        cls.getlinename(code).startswith(shellname)
                        for shellname in shellnames
                    )
                    for code in cls.decompose(linecode)
                )

            lines = [line for line in lines if valid(line.code)]

        # Energy selection
        if energybounds is not None:
            Z = energybounds[0]
            if not instance.isinteger(Z):
                Z = xraylib.SymbolToAtomicNumber(Z)

            def valid(energy):
                return (
                    energy >= energybounds[1]
                    and energy <= energybounds[2]
                    and energy != 0
                )

            lines = [line for line in lines if valid(line.energy(Z))]

        # Return
        return lines

    def __init__(self, line):
        """
        Args:
            code(int or str): xraylib line
        """
        if isinstance(line, self.__class__):
            self.code = line.code
            self.name = line.name
        else:
            self.code = self.getlinecode(line)
            self.name = self.getlinename(line)

    def __getstate__(self):
        return {"code": self.code, "name": self.name}

    def __setstate__(self, state):
        self.code = state["code"]
        self.name = state["name"]

    def _sortkey(self, other):
        return self.code

    def _cmpkey(self, other):
        if instance.isnumber(other):
            return self.code
        else:
            return repr(self)

    @property
    def _repr(self):
        """Unique representation of an instance"""
        return self.name

    def radrate(self, Z):
        """Radiative rate of a line: probability of this line / probabilty of fluorescence

        Args:
            Z(num): atomic number

        Returns:
            num
        """
        rate = xraylib.RadRate(Z, self.code)

        if rate == 0:
            # Maybe the radiative rate of one of its composites is known
            if self.code in xraylib.rcomposites:
                for comp in xraylib.rcomposites[self.code]:
                    if comp >= 0:
                        continue
                    rate = xraylib.RadRate(Z, comp)
                    if rate != 0:
                        # In abscence of a better assumption, assume all lines
                        # in the composite are equally probable
                        return rate / len(xraylib.composites[comp])

        return rate

    def energy(self, Z):
        """
        Args:
            Z(num): atomic number
        Returns:
            num
        """
        energy = xraylib.LineEnergy(Z, self.code)

        if energy == 0:
            # Maybe the energy of one of its composites is known
            if self.code in xraylib.rcomposites:
                for comp in xraylib.rcomposites[self.code]:
                    if comp >= 0:
                        continue
                    energy = xraylib.LineEnergy(Z, comp)
                    if energy != 0:
                        return energy

        return energy

    @property
    def shell(self):
        """Shell to which this line belongs to"""
        for shellname in xraylib.shell_to_code:
            if self.name.startswith(shellname):
                return Shell(shellname)
        raise RuntimeError("Cannot find the shell of fluorescence line {}".format(self))

    @property
    def groupname(self):
        """Group to which this line belongs to. Fluorescence lines are grouped by
        excitation shell, except the K-shell which is split in KA and KB
        """
        s = self.shell
        if self.code in xraylib.rcomposites:
            name = xraylib.code_to_line[xraylib.rcomposites[self.code][0]][0]
            if name == "KA" or name == "KB":
                return name
        return str(s)

    @property
    def shellsource(self):
        """Shell which is ionized after the transition"""
        for shellname in xraylib.shell_to_code:
            if self.name.endswith(shellname):
                return Shell(shellname)
        raise RuntimeError(
            "Cannot find the source shell of fluorescence line {}".format(self)
        )

    def fluorescence_production_cs(self, Z, E):
        # Kissel without cascade and Coster–Kronig transitions:
        #   return self.shell.partial_photo(Z,E)*self.shell.fluoyield(Z)*self.radrate(Z)
        # Kissel without cascade but with Coster–Kronig transitions:
        #   return xraylib.CS_FluorLine_Kissel_no_Cascade(Z, self.code, E)
        # Kissel with cascade and Coster–Kronig transitions:
        return xraylib.CS_FluorLine_Kissel_Cascade(Z, self.code, E)


class Shell(CompHashable, Copyable):
    @staticmethod
    def getshellname(shell):
        if isinstance(shell, Shell):
            return shell.name
        elif isinstance(shell, numbers.Number):
            return xraylib.code_to_shell[shell]
        else:
            if shell not in xraylib.shell_to_code.keys():
                raise ValueError("Unknown shell name {}".format(shell))
            return shell

    @staticmethod
    def getshellcode(shell):
        if isinstance(shell, Shell):
            return shell.code
        elif isinstance(shell, numbers.Number):
            if shell not in xraylib.code_to_shell.keys():
                raise ValueError("Unknown shell code {}".format(shell))
            return shell
        else:
            return xraylib.shell_to_code[shell]

    @classmethod
    def all_shells(cls, fluolines=None):
        shells = range(xraylib.shell_min, xraylib.shell_max + 1)
        return [cls(shell, fluolines=fluolines) for shell in shells]

    @classmethod
    def factory(cls, energybounds=None):
        if energybounds is None:
            return cls.all_shells()
        else:
            Z = energybounds[0]
            if not instance.isinteger(Z):
                Z = xraylib.SymbolToAtomicNumber(Z)

            def valid(energy):
                return energy >= energybounds[1] and energy <= energybounds[2]

            shells = cls.all_shells()
            shells = [s for s in shells if valid(s.edgeenergy(Z))]
            for s in shells:
                s._fluolines = FluoLine.factory(shells=[s], energybounds=energybounds)
            shells = [s for s in shells if s._fluolines]
            return shells

    @classmethod
    def expand(cls, shell):
        if isinstance(shell, Shell):
            return shell
        if instance.isnumber(shell):
            return shell
        if shell[1:]:
            if shell[1:].isdigit():
                return shell
            else:
                return shell[0]
        else:
            shell = shell[0].upper()
            n = ord(shell) - ord("K") + 1
            if n == 1:
                return shell
            else:
                return ["{}{}".format(shell, i) for i in range(1, 2 * n)]

    @classmethod
    def expandall(cls, shells):
        if not instance.isarray(shells):
            shells = [shells]
        return sorted(set(listtools.flatten(cls.expand(s) for s in shells)))

    @classmethod
    def expandselection(cls, shells, selection):
        if not instance.isarray(shells):
            shells = [shells]
        if not instance.isarray(selection):
            selection = [selection]
        for shell in selection:
            shells = [s for s in shells if not s != shell]
            shells += cls.expand(shell)
        return sorted(shells)

    @classmethod
    def shrink(cls, shells):
        shellgrps = set(str(s)[0].upper() for s in shells)
        for shell in shellgrps:
            subshells = cls.expand(shell)
            if all(s in shells for s in subshells):
                shells = [s for s in shells if not str(s).startswith(shell)]
                shells.append(shell)
        return sorted(shells)

    @classmethod
    def pymcafactory(cls, energybounds=None, split=None, allowsplitM=False):
        shells = cls.factory(energybounds=energybounds)
        shells = cls.shrink([str(s) for s in shells])

        # Only K,L and M
        lst = ["K", "L", "M"]
        shells = [s for s in shells if s[0] in lst]

        # Split
        if split:
            shells = cls.expandselection(shells, split)
        if not allowsplitM:
            shells = [s for s in shells if s[0] != "M" or s == "M"]

        return shells

    def __init__(self, shell, fluolines=None):
        """
        Args:
            code(int or str or Shell): shell code or name
            fluolines(Optional(array(int or str or Fluoline))): emission lines (None: all, []: explicit all)
        """
        self.code = self.getshellcode(shell)
        self.name = self.getshellname(shell)
        if isinstance(shell, Shell) and fluolines is None:
            self._fluolines = shell._fluolines
        else:
            if fluolines is None:
                self._fluolines = None  # all lines
            else:
                self._fluolines = FluoLine.factory(shells=[self], fluolines=fluolines)

    def __getstate__(self):
        return {"code": self.code, "name": self.name, "fluolines": self._fluolines}

    def __setstate__(self, state):
        self.code = state["code"]
        self.name = state["name"]
        self._fluolines = state["fluolines"]

    @property
    def fluolines(self):
        if self._fluolines is None:
            # None means all lines, but I need them explicitely now
            self._fluolines = FluoLine.factory(shells=[self])
        return self._fluolines

    def markinfo(self):
        if self._fluolines is None or self._fluolines == []:
            yield "{} fluorescence lines: all".format(self.name)
        else:
            yield "{} fluorescence lines: {}".format(
                self.name, ", ".join((str(line) for line in self._fluolines))
            )

    def _sortkey(self, other):
        return self.code

    def _cmpkey(self, other):
        if instance.isnumber(other):
            return self.code
        else:
            return repr(self)

    @property
    def _repr(self):
        """Unique representation of an instance"""
        return self.name

    def fluoyield(self, Z):
        """Fluorescence yield for this shell: probability for fluorescence / probability of shell ionization"""
        try:
            return xraylib.FluorYield(Z, self.code)
        except ValueError:
            return 0.0

    def radrate(self, Z):
        """Radiative rate of a shell: probabilities of selected lines / probabilty of fluorescence"""
        if self._fluolines is None:
            return [1]  # ALL lines
        else:
            return [line.radrate(Z) for line in self._fluolines]

    def atomiclevelwidth(self, Z):
        try:
            return xraylib.AtomicLevelWidth(Z, self.code)
        except ValueError:
            return 0.0

    def partial_fluoyield(self, Z, decomposed=False):
        """Probability of selected lines / probability of shell ionization"""
        if decomposed:
            fluoyield = self.fluoyield(Z)
            return {line: fluoyield * line.radrate(Z) for line in self.fluolines}
        else:
            if self._fluolines is None:
                return self.fluoyield(Z)
            else:
                return self.fluoyield(Z) * sum(self.radrate(Z))

    def partial_photo(self, Z, E):
        try:
            return xraylib.CS_Photo_Partial(Z, self.code, np.float64(E))
        except ValueError:
            return E * 0.0

    def edgeenergy(self, Z):
        return xraylib.EdgeEnergy(Z, self.code)


class FluoZGroup(CompHashable, Copyable):
    def __init__(self, element, group):
        self.element = element
        self.group = group
        self.element.markabsorber(shells=group)

    def _sortkey(self, other):
        return max(s.edgeenergy(self.element.Z) for s in self.element.shells)

    @property
    def _repr(self):
        """Unique representation of an instance"""
        return "{}-{}".format(self.element, self.group)


class FluoZLine(CompHashable, Copyable):
    def __init__(self, element, line):
        self.line = line
        self.element = element

    def __getattr__(self, attr):
        try:
            return getattr(self.line, attr)
        except Exception:
            return getattr(self.element, attr)

    @property
    def linewidth(self):
        """Lorentzian FWHM"""
        # Journal of Physical and Chemical Reference Data 8, 329 (1979); https://doi.org/10.1063/1.555595
        return self.line.shell.atomiclevelwidth(
            self.element.Z
        ) + self.line.shellsource.atomiclevelwidth(self.element.Z)

    def _sortkey(self, other):
        return self.energy()

    @property
    def _repr(self):
        """Unique representation of an instance"""
        return "{}-{}".format(self.element, self.line)

    @property
    def groupname(self):
        return "{}-{}".format(self.element, self.line.groupname)

    @property
    def nenergy(self):
        return 1

    def energy(self, **kwargs):
        return self.line.energy(self.element.Z)

    def split(self):
        return [self]

    @property
    def valid(self):
        return self.energy() > 0


class Line(CompHashable, Copyable):
    def __init__(self, energy, groupname="unknown", linewidth=0):
        self._energy = energy
        self._groupname = groupname
        self._linewidth = linewidth

    @property
    def linewidth(self):
        return self._linewidth

    def _sortkey(self, other):
        en = self._energy
        if instance.isiterable(en):
            en = max(en)
        return en

    @property
    def _repr(self):
        """Unique representation of an instance"""
        return str(self._energy)

    @property
    def nenergy(self):
        return listtools.length(self._energy)

    def energy(self, **kwargs):
        return self._energy

    @property
    def groupname(self):
        return self._groupname

    @property
    def valid(self):
        return True


class ScatteringLine(CompHashable, Copyable):
    def __init__(self, energysource):
        self.energysource = energysource

    @property
    def linewidth(self):
        return 0  # TODO: depends on many things

    def _sortkey(self, other):
        en = self.energy(polar=1e-3)
        if instance.isiterable(en):
            en = max(en)
        return en

    @property
    def _repr(self):
        """Unique representation of an instance"""
        return self.name

    @property
    def nenergy(self):
        return listtools.length(self.energysource)

    def split(self):
        if self.nenergy == 1:
            return [self]
        else:
            return [self.__class__(en) for en in self.energysource]

    @property
    def groupname(self):
        return str(self)

    @property
    def valid(self):
        return True


class RayleighLine(ScatteringLine):
    def __init__(self, energysource):
        self.name = "Rayleigh"
        super(RayleighLine, self).__init__(energysource)

    def energy(self, **kwargs):
        return self.energysource


class ComptonLine(ScatteringLine):
    def __init__(self, energysource):
        self.name = "Compton"
        super(ComptonLine, self).__init__(energysource)

    def energy(self, polar=None, **kwargs):
        """
        Args:
            polar(num): radians
        """
        if not polar:
            return self.energysource
        m = (
            ureg.Quantity(1 - np.cos(polar), "1/(m_e*c^2)")
            .to("1/keV", "spectroscopy")
            .magnitude
        )
        return self.energysource / (1 + self.energysource * m)

    def scatteringangle(self, energy):
        """
        Args:
            energy(num): scattering energy in keV

        Returns:
            polar(num): radians
        """
        m = (
            ureg.Quantity(1.0 / energy - 1.0 / self.energysource, "1/keV")
            .to("1/(m_e*c^2)", "spectroscopy")
            .magnitude
        )
        return np.arccos(np.clip(1 - m, -1, 1))


class Spectrum(Copyable, collections.abc.MutableMapping):

    TYPES = Enum(["diffcrosssection", "crosssection", "rate"])
    # diffcrosssection: cm^2/g/srad
    # crosssection: cm^2/g
    # rate: dimensionless

    def __init__(self, *args, **kwargs):
        self._lines = {}
        self.density = kwargs.pop("density", None)
        self.xlim = kwargs.pop("xlim", None)
        self.title = kwargs.pop("title", None)
        self.xlabel = "Energy (keV)"
        self.type = kwargs.pop("type", None)
        self.geometry = kwargs.pop("geometry", None)
        self.geomkwargs = kwargs.pop("geomkwargs", {})
        self.update(*args, **kwargs)

    def __getitem__(self, item):
        # Line values can be a function of geomkwargs
        value = self._lines[item]
        if callable(value):
            return value(**self.geomkwargs)
        else:
            return value

    def getitem_converted(self, item, convert=False, **kwargs):
        m = self.line_multiplier(convert=convert)
        return self[item] * m

    def __setitem__(self, item, value):
        self._lines[item] = value

    def __delitem__(self, item):
        del self._lines[item]

    def __iter__(self):
        return iter(self._lines)

    def __len__(self):
        return len(self._lines)

    @property
    def geomkwargs(self):
        if self.geometry is None:
            return self._geomkwargs
        else:
            return self.geometry.xrayspectrumkwargs()

    @geomkwargs.setter
    def geomkwargs(self, value):
        keys = ["polar", "azimuth"]
        if not hasattr(self, "_geomkwargs"):
            self._geomkwargs = {}
        for k in keys:
            self._geomkwargs[k] = value.get(k, None)

    @property
    def mcagain(self):
        if self.geometry is None:
            return 0.005
        else:
            return self.geometry.detector.mcagain

    @property
    def mcazero(self):
        if self.geometry is None:
            return 0
        else:
            return self.geometry.detector.mcazero

    def apply_weights(self, weights):
        if weights is None:
            return
        weights, func = instance.asarrayf(weights)
        n = len(weights)
        weights = func(weights / weights.sum(dtype=float))
        for k in self:
            if listtools.length(self[k]) == n:
                self[k] = self[k] * weights

    def sum_sources(self):
        for k in self:
            if isinstance(k, FluoZLine):
                self[k] = np.sum(self[k])

    def ylabel(self, convert=True, phsource=False, mcabin=False):
        if self.type == self.TYPES.crosssection:
            if convert:
                if phsource:
                    label = "Intensity"
                    units = "ph/cm/srad"
                else:
                    label = "Probability"
                    units = "1/cm/srad"
            else:
                label = "Cross-section"
                units = "cm$^2$/g"
                if phsource:
                    units = "{}.ph".format(units)
        elif self.type == self.TYPES.diffcrosssection:
            if convert:
                if phsource:
                    label = "Intensity"
                    units = "ph/cm/srad"
                else:
                    label = "Probability"
                    units = "1/cm/srad"
            else:
                label = "Differential cross-section"
                units = "cm$^2$/g/srad"
                if phsource:
                    units = "{}.phsource".format(units)
        else:
            if phsource:
                label = "Intensity"
                units = "ph"
            else:
                label = "Rate"
                units = "I/I$_0$"
        if mcabin:
            units = "{}/mcabin".format(units)
        return "{} {}".format(label, units)

    def line_multiplier(self, convert=True):
        if not convert:
            return 1
        if self.type == self.TYPES.crosssection:
            return self.density / (4 * np.pi)  # cm^2/g -> 1/cm/srad
        elif self.type == self.TYPES.diffcrosssection:
            return self.density  # cm^2/g/srad -> 1/cm/srad
        elif self.type == self.TYPES.rate:
            return 1
        else:
            raise RuntimeError("Unknown spectrum type: {}".format(self.type))

    def profile_multiplier(self, fluxtime=None, histogram=False):
        ret = 1
        if fluxtime is not None:
            ret = ret * fluxtime
        if histogram:
            ret = ret * self.mcagain
        return ret

    def __iadd__(self, other):
        if self.type != other.type:
            raise ValueError("Spectra need to be of the same type")
        for k, v in other.items():
            if k in self:
                self._lines[k] += v
            else:
                self._lines[k] = v
        return self

    def __add__(self, other):
        ret = copy.deepcopy(self)
        ret += other
        return ret

    def items_sorted(self, sort=False):
        if sort:

            def sortkey(x):
                return np.max(instance.asarray(x[0].energy(**self.geomkwargs)))

            return sorted(self.items(), key=sortkey)
        else:
            return self.items()

    def items_converted(self, convert=True, **kwargs):
        m = self.line_multiplier(convert=convert)
        for line, v in self.items_sorted(**kwargs):
            yield line, v * m

    @property
    def probabilities(self):
        """
        Interaction probability (1/cm/srad)
        """
        if (
            self.type != self.TYPES.crosssection
            and self.type != self.TYPES.diffcrosssection
        ):
            raise RuntimeError(
                "Spectrum does not contain cross-sections (cm^2/g or cm^2/g/srad)"
            )
        return self.items_converted(convert=True)

    @property
    def probabilities_dict(self):
        """
        Interaction probability (1/cm/srad)
        """
        ret = {}
        for k, v in self.probabilities:
            if k in ret:
                ret[k] += v
            else:
                ret[k] = v
        return ret

    @property
    def rates(self):
        """
        Line rate (ph/phsource)
        """
        if self.type != self.TYPES.rate:
            raise RuntimeError("Spectrum does not contain line rates (ph/phsource)")
        return self.items_converted(convert=True)

    @property
    def rates_dict(self):
        """
        Line rate (ph/phsource)
        """
        ret = {}
        for k, v in self.rates:
            if k in ret:
                ret[k] += v
            else:
                ret[k] = v
        return ret

    def lines(self, **kwargs):
        """
        Line rates (ph/phsource) or probabilities (1/cm/srad)
        """
        for line, v in self.items_converted(convert=True, **kwargs):
            if isinstance(line, FluoZLine):
                yield line, np.sum(v)
            else:
                energysource = instance.asarray(line.energy())
                values = instance.asarray(v)
                cls = line.__class__
                for e, v in zip(energysource, values):
                    yield cls(e), v

    def spectrum(self, **kwargs):
        """
        Line energy and rates (ph/phsource) or probabilities (1/cm/srad)
        """
        for line, v in self.lines(**kwargs):
            yield instance.asscalar(line.energy(**self.geomkwargs)), v

    @property
    def spectrum_dict(self):
        """
        Line energy and rates (ph/phsource) or probabilities (1/cm/srad)
        """
        return dict(Counter(self.spectrum()))

    def linegroups(self, **kwargs):
        ret = {}
        for line, v in self.items_converted(**kwargs):
            group = line.groupname
            if group not in ret:
                ret[group] = {}
            if isinstance(line, FluoZLine):
                # add contribution from all source lines
                ret[group][line] = np.sum(v)
            else:
                ret[group][line] = v
        return ret

    def lineinfo(self, convert=True, sort=True):
        lineinfo = {}
        geomkwargs = self.geomkwargs
        groups = self.linegroups(convert=convert)
        for group, lines in groups.items():
            bscat = group == "Compton" or group == "Rayleigh"
            if bscat:
                lines, intensities = next(iter(lines.items()))
                intensities = instance.asarray(intensities)
                lines = lines.split()
            else:
                intensities = list(lines.values())
            imax = np.argmax(intensities)
            for i, (line, peakarea) in enumerate(zip(lines, intensities)):
                energy = line.energy(**geomkwargs)
                if i == imax:
                    label = group
                else:
                    label = None
                if bscat:
                    linename = "{}{}".format(group, i)
                else:
                    linename = str(line)
                lineinfo[linename] = {
                    "group": group,
                    "label": label,
                    "energy": energy,
                    "area": peakarea,
                    "natwidth": line.linewidth,
                }
        if sort:
            lineinfo = collections.OrderedDict(
                sorted(lineinfo.items(), key=lambda x: x[1]["energy"])
            )
        return lineinfo

    @staticmethod
    def _lineinfo_values(lineinfo, key):
        return listtools.numpy_flatten(v[key] for v in lineinfo.values())

    def peakprofiles(self, x, lineinfo, voigt=False, **kwargs):
        if self.geometry is None:
            return None
        energies = self._lineinfo_values(lineinfo, "energy")
        if voigt:
            linewidths = self._lineinfo_values(lineinfo, "natwidth")
        else:
            linewidths = np.zeros_like(energies)
        return self.geometry.detector.lineprofile(
            x, energies, linewidth=linewidths, **kwargs
        )

    def peakprofiles_selectioninfo(self, lineinfo, lines, voigt=False):
        if self.geometry is None:
            raise RuntimeError("Cannot calculate detection limits without a geometry")

        if not instance.isarray(lines):
            lines = [lines]

        energies = self._lineinfo_values(lineinfo, "energy")
        areas = self._lineinfo_values(lineinfo, "area")
        if voigt:
            linewidths = self._lineinfo_values(lineinfo, "natwidth")
        else:
            linewidths = np.zeros_like(energies)

        # Indices of peaks in/out
        lineinfok = list(lineinfo.keys())
        indin = [lineinfok.index(k) for k in lines]
        indout = [k for k in range(len(lineinfo)) if k not in indin]
        info = {"indin": indin, "indout": indout}

        # Profile generators
        energiesin = energies[indin]
        linewidthsin = linewidths[indin]
        areasin = areas[indin]

        areasin = np.atleast_1d(areasin)[np.newaxis, :]
        areas = areas[np.newaxis, :]

        info["profilesin"] = lambda x: self.geometry.detector.lineprofile(
            x, energiesin, linewidth=linewidthsin, normalized=False, decomposed=True
        )
        info["profilesall"] = lambda x: self.geometry.detector.lineprofile(
            x, energies, linewidth=linewidths, normalized=False, decomposed=True
        )
        info["profiles"] = (
            lambda x: self.geometry.detector.lineprofile(
                x, energies, linewidth=linewidths, normalized=False, decomposed=False
            )
            * areas
        )

        # Extract from profiles:
        info["extractallpeaks"] = lambda profs: (profs[0] * areas).sum(axis=-1)
        info["extractinpeaks"] = lambda profs: (profs[0] * areasin).sum(axis=-1)

        # info["extractinpeaks_fromall"] = lambda profs: np.atleast_2d((profs[0]*areas)[:,indin]).sum(axis=-1)
        info["extractbkg"] = lambda profs: np.atleast_2d(
            (profs[0] * areas)[:, indout]
        ).sum(axis=-1) + (sum(profs[1:]) * areas).sum(axis=-1)

        # ROI generator
        std = np.sqrt(self.geometry.detector.voigtVAR(energiesin, linewidthsin))
        roi = zip(energiesin, std)

        def roigen(kstd=3, roi=roi):
            for r in roi:
                w = kstd * r[1]
                a, b = r[0] - w, r[0] + w
                yield a, b

        info["roigen"] = roigen

        return info

    def scaleprofiles(
        self, convert=False, fluxtime=None, histogram=False, lineinfo=None
    ):
        multiplier = self.profile_multiplier(fluxtime=fluxtime, histogram=histogram)
        ylabel = self.ylabel(
            convert=convert, phsource=fluxtime is not None, mcabin=histogram
        )
        if lineinfo is not None:
            for k in lineinfo:
                lineinfo[k]["area"] *= multiplier
        return multiplier, ylabel

    def __str__(self):
        lines = "\n ".join(
            ["{} {}".format(line, v) for line, v in self.lines(sort=True)]
        )
        return "{}\n Line   {}\n {}".format(
            self.title, self.ylabel(convert=True), lines
        )

    @property
    def energysource(self):
        return self["Rayleigh"].energy()

    @property
    def nsource(self):
        return listtools.length(self.energysource)

    def power(self, flux, **kwargs):
        if self.type != self.TYPES.rate:
            raise RuntimeError("Spectrum must contain rates, not cross-sections.")
        P = sum([energy * rate for energy, rate in self.spectrum(**kwargs)])
        return ureg.Quantity(P * flux, "keV/s")

    @property
    def channellimits(self):
        return self.energytochannel(self.xlim)

    def energytochannel(self, energy):
        energy, func = instance.asarrayf(energy)
        return func(
            np.clip(
                np.round((energy - self.mcazero) / self.mcagain).astype(int), 0, None
            )
        )

    @property
    def energylimits(self):
        return self.mcazero + self.mcagain * self.channellimits

    def mcachannels(self):
        a, b = self.channellimits
        return np.arange(a, b + 1)

    def mcaenergies(self):
        return self.mcazero + self.mcagain * self.mcachannels()

    def _linespectra(
        self, convert=True, fluxtime=None, histogram=False, voigt=False, sort=True
    ):
        energies = self.mcaenergies()
        lineinfo = self.lineinfo(convert=convert, sort=sort)

        # Normalized profiles: nchannels x npeaks
        profiles = self.peakprofiles(energies, lineinfo, voigt=voigt)
        if profiles is None:
            raise RuntimeError("Cannot calculate peak profiles without a geometry")

        # Real profiles: incoming photons and histogram
        _, ylabel = self.scaleprofiles(
            convert=convert, fluxtime=fluxtime, histogram=histogram, lineinfo=lineinfo
        )
        areas = self._lineinfo_values(lineinfo, "area")
        profiles *= areas[np.newaxis, :]

        return energies, profiles, ylabel, lineinfo

    def snr(
        self,
        linevalid,
        convert=True,
        fluxtime=None,
        histogram=False,
        voigt=False,
        backfunc=None,
        plot=False,
        kstd=3,
    ):
        energies = self.mcaenergies()
        lineinfo = self.lineinfo(convert=convert)
        lines = [k for k in lineinfo if linevalid(k)]
        if not bool(lines):
            raise RuntimeError("No lines fit the description")

        # Profile functions
        _, ylabel = self.scaleprofiles(
            convert=convert, fluxtime=fluxtime, histogram=histogram, lineinfo=lineinfo
        )
        info = self.peakprofiles_selectioninfo(lineinfo, lines, voigt=voigt)
        if backfunc is None:

            def backfunc(x):
                return np.zeros_like(x)

        def funcint(x):
            return info["extractinpeaks"](info["profilesin"](x))

        def functot(x):
            return np.squeeze(info["profiles"](x).sum(axis=-1) + backfunc(x))

        def funcbkg(x):
            return np.squeeze(info["extractbkg"](info["profilesall"](x)) + backfunc(x))

        # Plot spectrum
        if plot:
            lines = plt.plot(
                energies, np.random.poisson(np.round(functot(energies)).astype(int))
            )
            color2 = lines[0].get_color()
            color1 = next(plt.gca()._get_lines.prop_cycler)["color"]
            ax = plt.gca()

        # Signal and background
        total = 0.0
        background = 0.0
        arrtotal = []
        arrbackground = []
        for a, b in mergeroi1d(info["roigen"](kstd=kstd)):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                total += scipy.integrate.quad(functot, a, b)[0]
                background += scipy.integrate.quad(funcbkg, a, b)[0]
            if plot:
                ax.axvline(x=a)
                ax.axvline(x=b)
            a = np.argmin(np.abs(energies - a))
            b = np.argmin(np.abs(energies - b))
            ptotal = functot(energies[a:b])
            pbackground = funcbkg(energies[a:b])
            arrtotal.append(ptotal)
            arrbackground.append(pbackground)
            if plot:
                plt.fill_between(energies[a:b], pbackground, 0, alpha=0.5, color=color1)
                plt.fill_between(
                    energies[a:b], ptotal, pbackground, alpha=0.5, color=color2
                )
        signal = total - background
        noise = np.sqrt(total + background)
        SNR = signal / noise

        # Errors from linear fit within ROI:
        #   ytotal = a*signal + b*background
        arrtotal = np.concatenate(arrtotal)
        arrbackground = np.concatenate(arrbackground)
        arrsignal = arrtotal - arrbackground
        A = np.stack([arrsignal / signal, arrbackground / background], axis=-1)
        esignal, _ = fit1d.lstsq_std(A, vare=arrtotal)
        ESRdep = esignal / signal
        esignal, _ = fit1d.lstsq_std_indep(A, vare=arrtotal)
        ESRindep = esignal / signal
        cov = fit1d.lstsq_cov(A, vare=arrtotal)
        cor = fit1d.cor_from_cov(cov)
        cor = abs(cor[0, 1])
        if plot:
            plt.ylabel(ylabel)
            plt.xlabel("Energy (keV)")
            plt.title(
                "SNR = {:.02f}, ESR = {:.02f}% (dep), ESR = {:.02f}% (indep), COR = {:.02f}".format(
                    SNR, ESRdep * 100, ESRindep * 100, cor
                )
            )
        return {
            "SNR": SNR,
            "ESR (dependent)": ESRdep,
            "ESR (independent)": ESRindep,
            "SB-correlation": cor,
            "signal": signal,
            "noise": noise,
            "background": background,
        }

    def linespectra(self, sort=True, **kwargs):
        """X-ray spectrum decomposed in individual lines"""
        return self._linespectra(sort=sort, **kwargs)

    def groupspectra(self, sort=True, **kwargs):
        """X-ray spectrum decomposed in element-shell groups"""
        energies, profiles, ylabel, lineinfo = self._linespectra(sort=False, **kwargs)

        # Sort groups on maximum peak intensity
        groups = {}
        for line in lineinfo.values():
            if line["label"]:
                groups[line["group"]] = line["energy"]
        groups = zip(*sorted(groups.items(), key=lambda x: x[1]))[0]

        # Add lines per group
        lineinfo2 = collections.OrderedDict(
            (g, collections.OrderedDict()) for g in groups
        )

        nchan, npeaks = profiles.shape
        ret = np.zeros((nchan, len(groups)), dtype=profiles.dtype)
        for i, line in enumerate(lineinfo):
            ln = lineinfo[line]
            g = ln["group"]
            lineinfo2[g][line] = ln
            j = groups.index(g)
            ret[:, j] += profiles[:, i]

        return energies, ret, ylabel, lineinfo2

    def sumspectrum(self, backfunc=None, **kwargs):
        """Total X-ray spectrum"""
        energies, profiles, ylabel, lineinfo = self._linespectra(sort=False, **kwargs)
        profile = profiles.sum(axis=-1)
        if backfunc:
            profile = profile + backfunc(energies)
        return energies, profile, ylabel

    def plot(
        self,
        convert=False,
        fluxtime=None,
        mark=True,
        ylog=False,
        decompose=True,
        histogram=False,
        backfunc=None,
        voigt=False,
        forcelines=False,
        legend=True,
        marker=None,
        sumlabel="sum",
        title="",
    ):
        """X-ray spectrum or cross-section lines"""
        ax = plt.gca()
        if decompose:
            energies = self.mcaenergies()
            lines = self.lineinfo(convert=convert, sort=True)
            multiplier, ylabel = self.scaleprofiles(
                convert=convert, fluxtime=fluxtime, histogram=histogram
            )
            profiles = self.peakprofiles(
                energies, lines, voigt=voigt, onlyheight=forcelines
            )
            colors = {}
            for lineinfo in lines.values():
                g = lineinfo["group"]
                if g not in colors:
                    colors[g] = next(ax._get_lines.prop_cycler)["color"]
            blines = profiles is None or forcelines
            calcbkg = not (blines or backfunc is None)
            if calcbkg:
                bkg = backfunc(energies)
            else:
                bkg = 0
            sumprof = None
            off = 0
            for ind, (linename, lineinfo) in enumerate(lines.items()):
                color = colors[lineinfo["group"]]
                if profiles is None:
                    prof = lineinfo["area"] * multiplier
                else:
                    prof = lineinfo["area"] * multiplier * profiles[:, ind]
                    if backfunc is not None and forcelines:
                        off = backfunc(lineinfo["energy"])
                if legend:
                    label = lineinfo["label"]
                else:
                    label = None
                h = self._plotline(
                    lineinfo, energies, prof + bkg, color=color, off=off, label=label
                )
                if mark and lineinfo["label"] is not None:
                    plt.annotate(
                        linename,
                        xy=(lineinfo["energy"], h + off),
                        color=color,
                        clip_on=True,
                    )
                if not blines:
                    if sumprof is None:
                        sumprof = prof
                    else:
                        sumprof += prof
            if calcbkg:
                if legend:
                    label = sumlabel + "bkg"
                else:
                    label = None
                self._plotline(None, energies, bkg, label=label)
        else:
            energies, sumprof, ylabel = self.sumspectrum(
                convert=convert, fluxtime=fluxtime, histogram=histogram
            )
            if backfunc is None:
                bkg = 0
            else:
                bkg = backfunc(energies)
        if sumprof is not None:
            if legend:
                label = sumlabel
            else:
                label = None
            plt.plot(energies, sumprof + bkg, label=label, marker=marker)
        if self.geometry is not None and ylog:
            ax.set_yscale("log", base=10)
        if legend:
            plt.legend(loc="best")
        ax.set_xlabel(self.xlabel)
        ax.set_ylabel(ylabel)
        ax.set_xlim(*self.xlim)
        try:
            if title is not None:
                if not bool(title):
                    title = self.title
                plt.title(title)
        except UnicodeDecodeError:
            pass

    def _plotline(self, line, energies, prof, off=0, **kwargs):
        if listtools.length(prof) == 1:
            energy = line["energy"]
            h = instance.asscalar(prof)
            off = instance.asscalar(off)
            plt.plot([energy, energy], [off, h + off], **kwargs)
        else:
            y = prof
            h = max(y)
            plt.plot(energies, y, **kwargs)
        return h
