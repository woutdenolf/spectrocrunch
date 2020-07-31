# -*- coding: utf-8 -*-

from . import basetask
from ..utils import units


class Task(basetask.Task):
    """Task with an 'xrfgeometry' dependency
    """

    def _parameters_defaults(self):
        super(Task, self)._parameters_defaults()
        self.optional_parameters |= {
            "diodeI0gain",
            "diodeItgain",
            "xrf_positions",
            "referenceflux",
            "referencetime",
            "defaultexpotime",
            "samplecovers",
            "transmissionfilters",
        }

    @property
    def qxrfgeometry(self):
        if not hasattr(self, "_qxrfgeometry"):
            self._qxrfgeometry = None
        if self._qxrfgeometry is None:
            nxprocess = self.find_dependency(method="xrfgeometry")
            if nxprocess:
                geometry = nxprocess.results["geometry"].read(parse=True)
                geometry.diodeI0.gain = self._qxrf_parameter("diodeI0gain")
                geometry.diodeIt.gain = self._qxrf_parameter("diodeItgain")
                geometry.xrf_positions = self._qxrf_parameter("xrf_positions")
                value = self._qxrf_parameter("referenceflux", optional=True)
                if value is not None:
                    geometry.reference = units.Quantity(value, "Hz")
                value = self._qxrf_parameter("referencetime", optional=True)
                if value is not None:
                    geometry.defaultreferencetime = value
                value = self._qxrf_parameter("defaultexpotime", optional=True)
                if value is not None:
                    geometry.defaultexpotime = value
                value = self._qxrf_parameter("samplecovers", optional=True)
                if value:
                    geometry.addsamplecovers(value)
                value = self._qxrf_parameter("transmissionfilters", optional=True)
                if value:
                    geometry.addtransmissionfilters(value)
                self._qxrfgeometry = geometry
        return self._qxrfgeometry

    def _qxrf_parameter(self, name, optional=False):
        value = self.parameters.get(name, None)
        if value is None and not optional:
            raise basetask.MissingParameter(name)
        return value
