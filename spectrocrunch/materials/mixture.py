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

from . import multielementbase
from . import compoundfromlist
from . import element
from . import types
from . import stoichiometry
from ..utils import instance

import numpy as np
import fisx


class Mixture(multielementbase.MultiElementBase):

    def __init__(self, compounds, frac, fractype=None, name=None):
        """
        Args:
            compounds(list[obj]): list of compound objects
            frac(list[float]): compound fractions in the mixture
            fractype(types.fraction): compound fraction type
            name(Optional[str]): compound name
        """

        # Mixture name
        if name is None:
            self.name = str(id(self))
        else:
            self.name = name

        # Compound mole fractions
        fractions = np.asarray(frac)
        if fractype == types.fraction.mole:
            nfrac = fractions
        elif fractype == types.fraction.volume:
            MM = np.asarray([c.molarmass() for c in compounds])
            rho = np.asarray([c.density for c in compounds])
            nfrac = stoichiometry.frac_volume_to_mole(fractions, rho, MM)
        else:
            MM = np.asarray([c.molarmass() for c in compounds])
            nfrac = stoichiometry.frac_weight_to_mole(fractions, MM)

        # Compounds (no duplicates)
        self._compose_compounds(compounds, nfrac)

    def _compose_compounds(self, compounds, nfrac):
        self._compounds = {}
        for c, n in zip(compounds, nfrac):
            if c in self._compounds:
                self._compounds[c] += float(n)
            else:
                self._compounds[c] = float(n)

    def __getstate__(self):
        return {'name': self.name,
                'compounds_keys': list(self._compounds.keys()),
                'compounds_values': list(self._compounds.values())}

    def __setstate__(self, state):
        self.name = state['name']
        self._compounds = dict(zip(state['compounds_keys'],
                                   state['compounds_values']))

    @property
    def parts(self):
        return self._compounds

    def addcompounds(self, compounds, fractions, fractype):
        """Add compounds to the mixture

        Args:
            compounds(Compound or list): compounds
            fractions(num or list): compound fractions
            fractype(types.fraction): compound fraction type

        """
        if not instance.isarray(compounds):
            compounds = [compounds]
        if not instance.isarray(fractions):
            fractions = [fractions]
        fractions = np.asarray(fractions)
        if fractype == types.fraction.mole:
            nfrac = self.molefractions()
            compounds = list(nfrac.keys()) + compounds
            nfrac = np.asarray(list(nfrac.values()))
            nfrac = stoichiometry.add_frac(nfrac, fractions)
        elif fractype == types.fraction.volume:
            vfrac = self.volumefractions()
            compounds = list(vfrac.keys()) + compounds
            vfrac = np.asarray(list(vfrac.values()))
            vfrac = stoichiometry.add_frac(vfrac, fractions)
            MM = np.asarray([c.molarmass() for c in compounds])
            rho = np.asarray([c.density for c in compounds])
            nfrac = stoichiometry.frac_volume_to_mole(vfrac, rho, MM)
        else:
            wfrac = self.massfractions()
            compounds = list(wfrac.keys()) + compounds
            wfrac = np.asarray(list(wfrac.values()))
            wfrac = stoichiometry.add_frac(wfrac, fractions)
            MM = np.asarray([c.molarmass()for c in compounds])
            nfrac = stoichiometry.frac_weight_to_mole(wfrac, MM)
        self._compose_compounds(compounds, nfrac)

    def change_fractions(self, dfrac, fractype):
        """Change the compound fractions

        Args:
            frac(dict): compound fractions
            fractype(types.fraction): fraction type
        """
        compounds = list(dfrac.keys())
        fractions = np.asarray(list(dfrac.values()))
        if fractype == types.fraction.mole:
            nfrac = fractions
        elif fractype == types.fraction.volume:
            MM = np.asarray([c.molarmass() for c in compounds])
            rho = np.asarray([c.density for c in compounds])
            nfrac = stoichiometry.frac_volume_to_mole(fractions, rho, MM)
        else:
            MM = np.asarray([c.molarmass() for c in compounds])
            nfrac = stoichiometry.frac_weight_to_mole(fractions, MM)
        self._compose_compounds(compounds, nfrac)

    def tocompound(self, name=None):
        if not name:
            name = self.name
        tmp = self.elemental_molefractions()
        return compoundfromlist.CompoundFromList(list(tmp.keys()),
                                                 list(tmp.values()),
                                                 types.fraction.mole,
                                                 self.density,
                                                 name=name)

    def __str__(self):
        ws = self.massfractions()
        names, fractions = zip(*ws.items())
        names = list(map(str, names))
        names = ['(' + name + ')' if ' ' in name else name for name in names]
        return ' + '.join("{:.02f} wt% {}".format(w*100, name)
                          for name, w in zip(names, fractions))

    def __getitem__(self, compound):
        return self._compounds[compound]

    @property
    def elements(self):
        return list(set(e for c in self._compounds for e in c.elements))

    def molarmass(self, total=True):
        # equivalent to using molefractions when total=True
        # not equivalent to using molefractions when total=False
        nfrac = self.elemental_molefractions(total=total)
        MM = np.asarray([c.molarmass() for c in nfrac])
        nfrac = np.asarray(list(nfrac.values()))
        return (MM*nfrac).sum()

    def molarmasseff(self):
        return self.molarmass(total=False)

    @property
    def Zeff(self):
        # equivalent to using elemental_molefractions
        nfrac = self.molefractions(total=False)
        Z = np.asarray([c.Zeff for c in nfrac])
        nfrac = np.asarray(list(nfrac.values()))
        return (Z*nfrac).sum()

    def molefractions(self, total=True):
        if total:
            return dict(self._compounds)
        else:
            nfrac = np.asarray(list(self._compounds.values()))
            nfrac /= nfrac.sum()
            return dict(zip(self._compounds.keys(), nfrac))

    def massfractions(self):
        nfrac = self.molefractions()
        MM = np.asarray([c.molarmass() for c in nfrac])
        wfrac = stoichiometry.frac_mole_to_weight(
            np.asarray(list(nfrac.values())), MM)
        return dict(zip(nfrac.keys(), wfrac))

    def volumefractions(self):
        nfrac = self.molefractions()
        MM = np.asarray([c.molarmass() for c in nfrac])
        rho = np.asarray([c.density for c in nfrac])
        wfrac = stoichiometry.frac_mole_to_volume(
            np.asarray(list(nfrac.values())), rho, MM)
        return dict(zip(nfrac.keys(), wfrac))

    def fractions(self, fractype):
        if fractype == types.fraction.mole:
            return self.molefractions()
        elif fractype == types.fraction.volume:
            return self.volumefractions()
        else:
            return self.massfractions()

    @property
    def nelements(self):
        return sum(c.nelements for c in self._compounds)

    @property
    def density(self):
        nfrac = self.molefractions()
        MM = np.asarray([c.molarmass() for c in nfrac])
        rho = np.asarray([c.density for c in nfrac])
        nfrac = np.asarray(list(nfrac.values()))
        return stoichiometry.density_from_molefrac(nfrac, rho, MM)

    def elemental_molefractions(self, total=True):
        ret = {}
        c_nfrac = self.molefractions(total=total)
        for c in c_nfrac:
            e_nfrac = c.molefractions(total=total)
            for e in e_nfrac:
                if e in ret:
                    ret[e] += c_nfrac[c]*e_nfrac[e]
                else:
                    ret[e] = c_nfrac[c]*e_nfrac[e]
        return ret

    def elemental_massfractions(self):
        ret = {}
        c_wfrac = self.massfractions()
        for c in c_wfrac:
            e_wfrac = c.massfractions()
            for e in e_wfrac:
                if e in ret:
                    ret[e] += c_wfrac[c]*e_wfrac[e]
                else:
                    ret[e] = c_wfrac[c]*e_wfrac[e]
        return ret

    def _crosssection(self, method, E, fine=False, decomposed=False, **kwargs):
        """Calculate mixture cross-sections
        """
        bscat = self._cs_scattering(method)
        if bscat and not self.hasscatterers:
            if decomposed:
                return {}
            else:
                return E*0.
        c_wfrac = self.massfractions()
        if decomposed:
            ret = {}
            for c in c_wfrac:
                if bscat and not c.isscatterer:
                    continue
                ret[c] = {"w": c_wfrac[c], "cs": {}}
                if hasattr(c, 'structure') and fine:
                    environ = c
                else:
                    environ = None
                kwargs['environ'] = environ
                for e, w in c.massfractions().items():
                    cs = getattr(e, method)(E, **kwargs)
                    ret[c]["cs"][e] = {"w": w, "cs": cs}
        else:
            if self._cs_dict(method):
                ret = {}
                for c in c_wfrac:
                    if self._cs_scattering(method) and not c.isscatterer:
                        continue
                    if hasattr(c, 'structure') and fine:
                        environ = c
                    else:
                        environ = None
                    kwargs['environ'] = environ
                    e_wfrac = c.massfractions()
                    for e in c.parts:
                        cs = getattr(e, method)(E, **kwargs)
                        if not cs:
                            continue
                        w = c_wfrac[c]*e_wfrac[e]
                        for k, v in cs.items():
                            ret[k] = ret.get(k, 0) + w*v
            else:
                ret = None
                for c in c_wfrac:
                    if self._cs_scattering(method) and not c.isscatterer:
                        continue
                    if hasattr(c, 'structure') and fine:
                        environ = c
                    else:
                        environ = None
                    kwargs['environ'] = environ
                    e_wfrac = c.massfractions()
                    eret = sum(c_wfrac[c]*e_wfrac[e]*getattr(e, method)(E, **kwargs)
                               for e in c.parts)
                    if ret is None:
                        ret = eret
                    else:
                        ret += eret
        return ret

    def markinfo(self):
        yield "{}".format(self.name)
        for c in self._compounds:
            for s in c.markinfo():
                yield " {}".format(s)

    def markabsorber(self, symb=None, shells=None, fluolines=None, energybounds=None):
        """
        Args:
            symb(str): element symbol
        """
        for c in self._compounds:
            c.markabsorber(symb, shells=shells,
                           fluolines=fluolines, energybounds=energybounds)

    def unmarkabsorber(self):
        for c in self._compounds:
            c.unmarkabsorber()

    def hasabsorbers(self):
        return any([c.hasabsorbers() for c in self._compounds])

    def markscatterer(self, name=None):
        """
        Args:
            name(str): compound name
        """
        for c in self._compounds:
            c.markscatterer(name)

    def unmarkscatterer(self):
        for c in self._compounds:
            c.unmarkscatterer()

    def hasscatterers(self):
        return any([c.isscatterer() for c in self._compounds])

    def get_energy(self, energyrange, defaultinc=1):
        """Get absolute energies (keV) from a relative energy range (eV)
        """
        for c in self._compounds:
            ret = c.get_energy(energyrange, defaultinc=defaultinc)
            if ret is not None:
                return ret
        return None

    def topymca(self, cfg, defaultthickness=1e-4):
        r = self.massfractions()
        massfractions = list(r.values())
        names = [mat.topymca(cfg, defaultthickness=defaultthickness)
                 for mat in r]

        matname = self.pymcaname
        cfg["materials"][matname] = {'Comment': self.pymcacomment,
                                     'CompoundFraction': massfractions,
                                     'Thickness': defaultthickness,
                                     'Density': self.density,
                                     'CompoundList': names}
        return matname

    def tofisx(self, cfg, defaultthickness=1e-4):
        r = self.massfractions()
        massfractions = list(r.values())
        names = [mat.tofisx(cfg, defaultthickness=defaultthickness)
                 for mat in r]

        matname = self.pymcaname
        o = fisx.Material(matname, self.density,
                          defaultthickness, self.pymcacomment)
        o.setCompositionFromLists(names, massfractions)
        cfg.addMaterial(o, errorOnReplace=False)

        return matname
