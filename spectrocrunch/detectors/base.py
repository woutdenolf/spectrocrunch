# -*- coding: utf-8 -*-

from ..utils import instance
from ..utils import units
from ..math import noisepropagation
from ..materials import element
from ..utils.copyable import Copyable

import numpy as np
import fisx
import re


class Material(Copyable):

    DETMATERIALLABEL = "Detector"
    ATMOSPHERELABEL = "Atmosphere"

    def __init__(self, attenuators=None):
        # TODO: separate attenuators, beamfilters and detector material?
        if attenuators is None:
            attenuators = {}
        self.attenuators = attenuators

    def __getstate__(self):
        return {'attenuators': self._attenuators}

    def __setstate__(self, state):
        self._attenuators = state['attenuators']

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.attenuators == other.attenuators
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def addattenuator(self, name, material, thickness):
        thickness = units.umagnitude(thickness, 'cm')
        self._attenuators[name] = {"material": material,
                                   "thickness": thickness}

    def addsamplecover(self, name, material, thickness):
        self.addattenuator(name, material, thickness)
        self.addbeamfilter(None, material, thickness)

    def addbeamfilter(self, name, material, thickness):
        if not name:
            name = self.nextbeamfilter()
        elif not self.isbeamfilter(name):
            raise ValueError('This is not a beam filter name')
        self.addattenuator(name, material, thickness)

    def adddetectorfilter(self, name, material, thickness):
        if self.isbeamfilter(name):
            raise ValueError('This is a beam filter name')
        self.addattenuator(name, material, thickness)

    def removeattenuator(self, name):
        self._attenuators.pop(name, None)

    def removebeamfilters(self):
        names = [att for att in self.attenuators if self.isbeamfilter(att)]
        for name in names:
            self.removeattenuator(name)

    def removedetectorfilters(self):
        names = [att for att in self.attenuators if self.isdetectorfilter(att)]
        for name in names:
            self.removeattenuator(name)

    @property
    def attenuators(self):
        try:
            mat = self.geometry.atmosphere
            d = self.geometry.distance
            if mat and d:
                self.addattenuator(self.ATMOSPHERELABEL, mat,
                                   units.umagnitude(d, "cm"))
        except AttributeError:
            pass
        return self._attenuators

    @attenuators.setter
    def attenuators(self, value):
        self._attenuators = value

    def __str__(self):
        s1 = "\n ".join(["{} = {} cm".format(
            attinfo["material"], attinfo["thickness"]) for attinfo in self.beamfilters()])
        s2 = "\n ".join(["{} = {} cm".format(
            attinfo["material"], attinfo["thickness"]) for attinfo in self.detectorfilters()])
        if s1:
            s1 = "Beam filters:\n {}".format(s1)
        else:
            s1 = "Beam filters: None"
        if s2:
            s2 = "Detector filters:\n {}".format(s2)
        else:
            s2 = "Detector filters: None"
        if self.hasmaterial:
            s3 = " Material = {} ({} cm)".format(self.material, self.thickness)
        else:
            s3 = " Material = None"
        return "\n".join([s3, s1, s2])

    @property
    def hasmaterial(self):
        return self.DETMATERIALLABEL in self.attenuators

    @property
    def material(self):
        return self.attenuators[self.DETMATERIALLABEL]["material"]

    @property
    def thickness(self):
        return self.attenuators[self.DETMATERIALLABEL]["thickness"]

    @thickness.setter
    def thickness(self, value):
        self.attenuators[self.DETMATERIALLABEL]["thickness"] = value

    def addtopymca(self, setup, cfg):
        for attlabel, attinfo in self.attenuators.items():
            self.addtopymca_attenuator(
                setup, cfg, attlabel, attinfo["material"], attinfo["thickness"])
        attenuators = {}
        for attlabel, values in cfg["attenuators"].items():
            if values[0]:
                attenuators[attlabel] = values
        cfg["attenuators"] = attenuators

    def addtopymca_attenuator(self, setup, cfg, attlabel, material, thickness):
        matname = setup.addtopymca_material(
            cfg, material, defaultthickness=thickness)
        # Can only handle explicite single elements
        if attlabel == self.DETMATERIALLABEL:
            if isinstance(material, element.Element):
                matname = '{}1'.format(material)
        cfg["attenuators"][attlabel] = [
            1, matname, material.density, thickness, 1.0]

    def loadfrompymca(self, setup, cfg):
        self.attenuators = {}
        for attlabel, attinfo in cfg["attenuators"].items():
            if attlabel != "Matrix":
                self.loadfrompymca_attenuator(setup, cfg, attlabel, attinfo)

    def loadfrompymca_attenuator(self, setup, cfg, attlabel, attinfo):
        enabled, matname, density, thickness = attinfo[:4]
        if enabled:
            material = setup.loadfrompymca_material(cfg, matname, density)
            self.attenuators[attlabel] = {
                "material": material, "thickness": thickness}

    def isdetectorfilter(self, name):
        return not self.isbeamfilter(name) and name != self.DETMATERIALLABEL

    def beamfilters(self):
        return [self.attenuators[att] for att in self.attenuators
                if self.isbeamfilter(att)]

    def detectorfilters(self):
        return [self.attenuators[att] for att in self.attenuators
                if self.isdetectorfilter(att)]

    def match_beamfilter(self, name):
        return re.compile(r'^beamfilter(\d+)$').match(name.lower())

    def isbeamfilter(self, name):
        return bool(self.match_beamfilter(name))

    def nextbeamfilter(self):
        i = 0
        for att in self.attenuators:
            m = self.match_beamfilter(att)
            if m:
                i = max(i, int(m.groups()[0]))
        return 'BeamFilter{}'.format(i)

    def addtofisx(self, setup, cfg):
        for attlabel, attinfo in self.attenuators.items():
            cfg.addtofisx_material(attinfo["material"])

        atts = [[attinfo["material"].pymcaname, attinfo["material"].density,
                 attinfo["thickness"], 1.0] for attinfo in self.beamfilters()]
        if atts:
            setup.setBeamFilters(atts)

        atts = [[attinfo["material"].pymcaname, attinfo["material"].density,
                 attinfo["thickness"], 1.0] for attinfo in self.detectorfilters()]
        if atts:
            setup.setAttenuators(atts)

    def filter_absorbance(self, energy, source=False):
        """Absorbance by filters (before or after sample)
        """
        if source:
            atts = self.beamfilters()
        else:
            atts = self.detectorfilters()
        if atts:
            return sum([att["material"].mass_att_coeff(energy)*(att["thickness"]*att["material"].density) for att in atts])
        else:
            return np.zeros_like(energy, dtype=float)

    def filter_transmission(self, energy, source=False):
        """Transmission through filters (before or after sample)
        """
        return np.exp(-self.filter_absorbance(energy, source=source))

    def filter_attenuation(self, energy, source=False):
        """Attenuation by filters (before or after sample)
        """
        return 1-self.filter_transmission(energy, source=source)

    def absorbance(self, energy):
        """Absorbance by detector material
        """
        if self.hasmaterial:
            return self.material.mass_att_coeff(energy)*(self.thickness*self.material.density)
        else:
            return np.full_like(3, np.inf, dtype=float)

    def transmission(self, energy):
        """Transmission through detector material
        """
        return np.exp(-self.absorbance(energy))

    def attenuation(self, energy):
        """Attenuation by detector material
        """
        return 1-self.transmission(energy)

    def efficiency(self, energysource, energydet, withdetectorattenuation=True):
        """Filter+Detector attenuation

          T(energysource)*T(energydet)*A(energydet)

          T: transmission by filters (before sample and detector)
          A: attenuation by the detector material (withdetectorattenuation==True)

        Args:
            energysource: n0
            energydet: n1
            withdetectorattenuation(Optional(bool)): with detector attenuation or not

        Returns:
            array: n0 x n1
        """
        energysource = instance.asarray(energysource)
        energydet = instance.asarray(energydet)
        T0 = self.filter_transmission(energysource, source=True)
        T1 = self.filter_transmission(energydet, source=False)
        if withdetectorattenuation:
            A = self.attenuation(energydet)
        else:
            A = 1.
        return T0[:, np.newaxis]*(T1*A)[np.newaxis, :]


class SolidState(Material):

    def __init__(self, ehole=None, **kwargs):
        self.ehole = ehole
        super(SolidState, self).__init__(**kwargs)

    @property
    def ehole(self):
        return self._ehole

    @ehole.setter
    def ehole(self, value):
        if value is None:
            self._ehole = None
        else:
            self._ehole = units.Quantity(value, 'eV').to('eV')

    def __getstate__(self):
        state = super(SolidState, self).__getstate__()
        state['ehole'] = self.ehole
        return state

    def __setstate__(self, state):
        super(SolidState, self).__setstate__(state)
        self.ehole = state['ehole']

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            if not (super(SolidState, self).__eq__(other)):
                return False
            return self.ehole == other.ehole
        else:
            return False

    def __str__(self):
        return " Ionization energy = {:~}\n{}".format(self.ehole, super(SolidState, self).__str__())


class CentricCone(SolidState):

    def __init__(self, activearea=None, **kwargs):
        self.activearea = activearea  # cm^2
        super(CentricCone, self).__init__(**kwargs)

    def __getstate__(self):
        state = super(CentricCone, self).__getstate__()
        state['activearea'] = self.activearea
        return state

    def __setstate__(self, state):
        super(CentricCone, self).__setstate__(state)
        self.activearea = state['activearea']

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            if not (super(CentricCone, self).__eq__(other)):
                return False
            return self.activearea == other.activearea
        else:
            return False

    def __str__(self):
        return " Active area = {:~}\n{}".format(
            self.activearea, super(CentricCone, self).__str__())

    @property
    def activearea_rv(self):
        return self._activearea

    @property
    def activearea(self):
        return noisepropagation.E(self.activearea_rv)

    @activearea.setter
    def activearea(self, value):
        if value is None:
            self._activearea = None
        else:
            self._activearea = units.Quantity(value, 'cm^2').to('cm^2')

    @classmethod
    def solidangle_calc(cls, activearea=None, distance=None, solidangle=None):
        # Cone with source on-axis:
        #  solidangle = 2.pi.(1-cos(theta)) where theta the opening angle
        #  cos(theta) = d/sqrt(r^2+d^2) where d the source-detector and r the radius of the active area (assume disk)
        if solidangle is None:
            d2 = distance**2.
            r2 = activearea/np.pi
            costheta = distance/(r2+d2)**0.5
            return 2.*np.pi*(1.-costheta.to(units.dimensionless).magnitude)
        elif distance is None:
            r2 = activearea/np.pi
            c2 = (1-solidangle/(2.*np.pi))**2
            return (r2*c2/(1-c2))**0.5
        elif activearea is None:
            c2 = (1-solidangle/(2.*np.pi))**2
            return (1-c2)*distance**2*np.pi/c2
        else:
            raise RuntimeError(
                "Either distance, active area or solid angle must be None")

    def addtofisx(self, setup, cfg):
        super(CentricCone, self).addtofisx(setup, cfg)
        if self.hasmaterial:
            detector = fisx.Detector(
                self.material.name, self.material.density, self.thickness)
            detector.setActiveArea(self.activearea.to("cm^2").magnitude)
            detector.setDistance(self.geometry.distance.to("cm").magnitude)
            # detector.setMaximumNumberOfEscapePeaks(0)
            setup.setDetector(detector)

    def addtopymca(self, setup, cfg):
        super(CentricCone, self).addtopymca(setup, cfg)
        cfg["concentrations"]["area"] = self.activearea.to("cm^2").magnitude

    def loadfrompymca(self, setup, cfg):
        super(CentricCone, self).loadfrompymca(setup, cfg)
        self.activearea = cfg["concentrations"]["area"]
