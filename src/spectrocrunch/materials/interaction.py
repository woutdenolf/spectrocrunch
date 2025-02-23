from ..patch.pint import ureg
from ..utils.hashable import CompHashable
from . import element

import numpy as np


class Interaction(CompHashable):
    def __init__(self, name, energy, prob):
        self._name = name
        self._energy = energy
        self._prob = prob

    @property
    def energy(self):
        return self._energy

    def _sortkey(self, other):
        return self.energy

    @property
    def _repr(self):
        """Unique representation of an instance"""
        return self._name


class InteractionSource(Interaction):
    def __init__(self, energy, index):
        name = "Source-{}".format(index)
        prob = 1
        super(InteractionSource, self).__init__(name, energy, prob)


class InteractionFluo(Interaction):
    def __init__(self, el, shell, line):
        """
        Args:
            el(element or num or str):
            shell(Shell):
            line(FluoLine):
        """
        if not isinstance(el, element.Element):
            el = element.Element(el)
        name = "{}-{}".format(el, line)
        energy = line.energy(el.Z)
        prob = shell.fluoyield(el.Z) * line.radrate(el.Z)
        super(InteractionFluo, self).__init__(name, energy, prob)


class InteractionElScat(Interaction):
    def __init__(self, source):
        name = "RScat({})".format(source)
        prob = 1
        super(InteractionElScat, self).__init__(name, source.energy, prob)


class InteractionInelScat(Interaction):
    def __init__(self, source, theta):
        name = "CScat({})".format(source)
        prob = 1
        self.theta = theta  # scattering angle

        super(InteractionInelScat, self).__init__(name, source.energy, prob)

    @property
    def energy(self):
        if self.theta == 0:
            return self.energy
        delta = (
            ureg.Quantity(1 - np.cos(np.radians(self.theta)), "1/(m_e*c^2)")
            .to("1/keV", "spectroscopy")
            .magnitude
        )
        return self._energy / (1 + self._energy * delta)
