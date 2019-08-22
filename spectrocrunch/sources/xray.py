# -*- coding: utf-8 -*-

from ..utils.classfactory import with_metaclass
from ..utils.Enum import Enum
from . import polarization
from ..utils.copyable import Copyable

import numpy as np
import matplotlib.pyplot as plt


class XraySource(with_metaclass(Copyable)):

    def __init__(self, stokes=None):
        self.stokes = stokes

    def __getstate__(self):
        return {'stokes': self.stokes}

    def __setstate__(self, state):
        self.stokes = state['stokes']

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.stokes == other.stokes
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __str__(self):
        s = str(self.stokes).replace('\n', '\n ')
        return "XraySource:\n {}".format(s)

    @property
    def thomson_K(self):
        return self.stokes.thomson_K

    def compton_K(self, energy):
        return self.stokes.compton_K(energy)

    def addtopymca(self, setup, cfg):
        pass

    def loadfrompymca(self, setup, cfg):
        pass


class Synchrotron(XraySource):

    def __init__(self, **polparams):
        if "intensity" not in polparams:
            polparams["intensity"] = 1  # W/m^2
        if "dop" not in polparams:
            polparams["dop"] = 1.  # degree of polarization (in [0,1])
        if "dolp" not in polparams:
            # degree of linear polarization (in [0,dop])
            polparams["dolp"] = 1.*polparams["dop"]
        if "polangle" not in polparams:
            # angle of polarization ellipse with respect to the horizontal direction (in [-90,90])
            polparams["polangle"] = 0
        # if "handedness" not in polparams:
        #    polparams["handedness"] = "left" # above/below the plane
        stokes = polarization.Stokes.from_params(**polparams)
        super(Synchrotron, self).__init__(stokes=stokes)


class Tube(XraySource):

    def __init__(self, intensity=1):
        stokes = polarization.Stokes.from_params(intensity=intensity,
                                                 dop=0,
                                                 dolp=0,
                                                 polangle=0,
                                                 handedness="left")
        super(Tube, self).__init__(stokes=stokes)


factory = XraySource.factory
registry = XraySource.clsregistry
