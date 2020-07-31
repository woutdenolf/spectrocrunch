# -*- coding: utf-8 -*-

from ..utils import lut
from ..utils import units
from ..utils.copyable import Copyable


class Optics(Copyable):
    def __init__(self, uselut=True, default=1, **kwargs):
        if uselut:
            self.lut = lut.LUT(default=default, **kwargs)
        self.uselut = uselut

    def __getstate__(self):
        state = {"uselut": self.uselut}
        if self.uselut:
            state["lut"] = self.lut
        return state

    def __setstate__(self, state):
        self.uselut = state["uselut"]
        if "lut" in state:
            self.lut = state["lut"]

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            if self.uselut ^ other.uselut:
                return False
            if self.uselut:
                return self.lut == other.lut
            else:
                return True
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __str__(self):
        name = type(self).__name__
        s = "\n ".join(
            "{:~}: {} %".format(k, v * 100) for k, v in self.lut.zip("keV", None)
        )
        if s:
            return "{}:\n transmission:\n {}".format(name, s)
        else:
            return "{}:\n transmission: 100%".format(name)

    def reset_transmission(self):
        if self.uselut:
            self.lut.clear(1)

    def transmission(self, energy):
        self.checklut()
        energy = units.Quantity(energy, "keV")
        T = self.lut(energy)
        return T.to(units.dimensionless).magnitude

    def set_transmission(self, energy, transmission):
        self.checklut()
        self.lut.replace(units.Quantity(energy, "keV"), transmission)

    def checklut(self):
        if not self.uselut:
            raise RuntimeError(
                "{} has no transmission lookup table.".format(type(self).__name__)
            )
