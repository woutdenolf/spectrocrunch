# -*- coding: utf-8 -*-

from . import types
from . import elementbase
from ..utils import instance

import numpy as np


class MultiElementBase(elementbase.ElementBase):
    # TODO: refactor most of Mixture and Compound

    @property
    def _repr(self):
        return self.name

    @property
    def nparts(self):
        return len(self.parts)

    def setfraction(self, parts, values, fractype):
        if instance.isstring(parts):
            parts = [parts]
        values = instance.asarray(values)

        for p in parts:
            if p not in self.parts:
                raise RuntimeError("{} not in {}".format(p, self))

        # rebalance others
        w = self.fractions(fractype)
        w2 = dict(w)
        for p in parts:
            w2.pop(p)

        # update others
        if w2:
            v2 = np.asarray(w2.values())
            v2 *= (1-values.sum())/v2.sum()
            w.update((k, v) for k, v in zip(w2.keys(), v2))

        # update fractions
        w.update(zip(parts, values))
        self.change_fractions(w, fractype)

    def setmassfraction(self, comp, value):
        self.setfraction(comp, value, types.fraction.mass)

    def setmolefraction(self, comp, value):
        self.setfraction(comp, value, types.fraction.mole)

    def setvolumefraction(self, comp, value):
        self.setfraction(comp, value, types.fraction.volume)

    @classmethod
    def _cs_scattering(cls, method):
        return method == "scattering_cross_section" or method == "compton_cross_section" or method == "rayleigh_cross_section" or\
            method == "diff_compton_cross_section" or method == "diff_rayleigh_cross_section"

    @classmethod
    def _cs_dict(cls, method):
        return method == "fluorescence_cross_section_lines" or method == "diff_fluorescence_cross_section"

    @classmethod
    def _cs_lazy(cls, method):
        return method == "diff_compton_cross_section" or method == "diff_rayleigh_cross_section"

    def mass_att_coeff(self, E, fine=False, decomposed=False, **kwargs):
        """Mass attenuation coefficient (cm^2/g, E in keV). Use for transmission XAS.
        """
        return self._crosssection("mass_att_coeff", E, fine=fine, decomposed=decomposed, **kwargs)

    def mass_abs_coeff(self, E, fine=False, decomposed=False, **kwargs):
        """Mass absorption coefficient (cm^2/g, E in keV).
        """
        return self._crosssection("mass_abs_coeff", E, fine=fine, decomposed=decomposed, **kwargs)

    def partial_mass_abs_coeff(self, E, fine=False, decomposed=False, **kwargs):
        """Mass absorption coefficient for the selected shells and lines (cm^2/g, E in keV).
        """
        return self._crosssection("partial_mass_abs_coeff", E, fine=fine, decomposed=decomposed, **kwargs)

    def scattering_cross_section(self, E, fine=False, decomposed=False, **kwargs):
        """Scattering cross section (cm^2/g, E in keV).
        """
        return self._crosssection("scattering_cross_section", E, fine=fine, decomposed=decomposed, **kwargs)

    def compton_cross_section(self, E, fine=False, decomposed=False, **kwargs):
        """Compton cross section (cm^2/g, E in keV).
        """
        return self._crosssection("compton_cross_section", E, fine=fine, decomposed=decomposed, **kwargs)

    def rayleigh_cross_section(self, E, fine=False, decomposed=False, **kwargs):
        """Rayleigh cross section (cm^2/g, E in keV).
        """
        return self._crosssection("rayleigh_cross_section", E, fine=fine, decomposed=decomposed, **kwargs)

    def fluorescence_cross_section(self, E, fine=False, decomposed=False, **kwargs):
        """XRF cross section (cm^2/g, E in keV). Use for fluorescence XAS.
        """
        return self._crosssection("fluorescence_cross_section", E, fine=fine, decomposed=decomposed, **kwargs)

    def fluorescence_cross_section_lines(self, E, fine=False, decomposed=False, **kwargs):
        """XRF cross section (cm^2/g, E in keV). Use for XRF.
        """
        return self._crosssection("fluorescence_cross_section_lines", E, fine=fine, decomposed=decomposed, **kwargs)

    def diff_fluorescence_cross_section(self, E, fine=False, decomposed=False, **kwargs):
        """Differential XRF cross section (cm^2/g/srad, E in keV). Use for XRF.
        """
        return self._crosssection("diff_fluorescence_cross_section", E, fine=fine, decomposed=decomposed, **kwargs)

    def diff_rayleigh_cross_section(self, E, fine=False, decomposed=False, **kwargs):
        """Differential Rayleigh cross section (cm^2/g/srad, E in keV). Use for XRF.
        """
        return self._crosssection("diff_rayleigh_cross_section", E, fine=fine, decomposed=decomposed, **kwargs)

    def diff_compton_cross_section(self, E, fine=False, decomposed=False, **kwargs):
        """Differential Compton cross section (cm^2/g/srad, E in keV). Use for XRF.
        """
        return self._crosssection("diff_compton_cross_section", E, fine=fine, decomposed=decomposed, **kwargs)

    @classmethod
    def cs_type(cls, cs):
        # 0. cs = [...]                -> pure element cs
        # 1. cs = {'w':0.1,'cs':[...]} -> element w+cs
        # 2. cs = {'A':{},'B':{}}      -> compound or mixture
        if isinstance(cs, dict):
            if 'cs' in cs:
                if isinstance(cs['cs'], dict):
                    return 2
                else:
                    return 1
            else:
                return 2
        else:
            return 0

    @classmethod
    def cs_collapse(cls, cs):
        t = cls.cs_type(cs)
        if t == 0:
            return cs
        elif t == 1:
            return cs['w']*cs['cs']
        else:
            return sum(cls.cs_collapse(c) for c in cs.values())

    @classmethod
    def csdict_parse(cls, cs):
        t = cls.cs_type(cs)
        if t == 0:
            return np.array([1]), [cs]
        elif t == 1:
            return np.array([cs['w']]), [cs['cs']]
        else:
            w = np.array([c['w'] for c in cs.values()])
            cs = [cls.cs_collapse(c['cs']) for c in cs.values()]
            return w, cs

    @classmethod
    def csdict_parse_elements(cls, csin, csout=None, w=1, k=None):
        if csout is None:
            csout = {}
        t = cls.cs_type(csin)
        if t == 0:
            csout[k] = csout.get(k, 0)+w*cs
        elif t == 1:
            csout[k] = csout.get(k, 0)+w*cs['w']*cs['cs']
        else:
            for k, c in csin.items():
                cls.csdict_parse_elements(
                    c['cs'], csout=csout, w=w*c['w'], k=k)
        return csout
