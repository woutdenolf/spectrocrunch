import numpy as np
from ..utils import units
from ..utils.copyable import Copyable
from ..math import noisepropagation


class Base(Copyable):
    def __init__(self, detector=None, source=None, atmosphere=None):
        self.detector = detector
        self.source = source
        self.atmosphere = atmosphere

    def __getstate__(self):
        return {
            "detector": self.detector,
            "source": self.source,
            "atmosphere": self.atmosphere,
        }

    def __setstate__(self, state):
        self.detector = state["detector"]
        self.source = state["source"]
        self.atmosphere = state["atmosphere"]

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            if not (self.detector == other.detector):
                return False
            if not (self.source == other.source):
                return False
            return self.atmosphere == other.atmosphere
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    @property
    def detector(self):
        return self._detector

    @detector.setter
    def detector(self, value):
        self._detector = value
        if self._detector is not None:
            self._detector.geometry = self

    def __getattr__(self, attr):
        try:
            return getattr(self._detector, attr)
        except AttributeError:
            try:
                return getattr(self.source, attr)
            except AttributeError:
                raise AttributeError(
                    "'{}' object has no attribute '{}'".format(
                        self.__class__.__name__, attr
                    )
                )

    def __str__(self):
        return "{}\n{}".format(self.source, self.detector)

    def addtofisx(self, setup, cfg):
        self.detector.addtofisx(setup, cfg)

    def addtopymca(self, setup, cfg):
        self.detector.addtopymca(setup, cfg)
        self.source.addtopymca(setup, cfg)

    def loadfrompymca(self, setup, cfg):
        self.detector.loadfrompymca(setup, cfg)
        self.source.loadfrompymca(setup, cfg)


class FlatSample(Base):
    def __init__(
        self,
        anglein=None,
        angleout=None,
        sample_azimuth=180.0,
        detector_azimuth=0.0,
        **kwargs,
    ):
        """
        Args:
            anglein(num): angle (deg) between primary beam and sample surface [0, 90]
            angleout(num): angle (deg) between detector and sample surface [-90, 90]  (>0 means reflection)
            sample_azimuth(num): angle (deg) for surface normal (pointing inwards) with the polarization plane
            detector_azimuth(num): angle (deg) for detector with the polarization plane
        """
        self.anglein = anglein  # deg
        self.angleout = angleout  # deg
        self.sample_azimuth = sample_azimuth  # deg
        self.detector_azimuth = detector_azimuth  # deg
        super(FlatSample, self).__init__(**kwargs)

    def __getstate__(self):
        state = super(FlatSample, self).__getstate__()
        state["anglein"] = self.anglein
        state["angleout"] = self.angleout
        state["sample_azimuth"] = self.sample_azimuth
        state["detector_azimuth"] = self.detector_azimuth
        return state

    def __setstate__(self, state):
        super(FlatSample, self).__setstate__(state)
        self.anglein = state["anglein"]
        self.angleout = state["angleout"]
        self.sample_azimuth = state["sample_azimuth"]
        self.detector_azimuth = state["detector_azimuth"]

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            if not (super(FlatSample, self).__eq__(other)):
                return False
            return (
                self.anglein == other.anglein
                and self.angleout == other.angleout
                and self.sample_azimuth == other.sample_azimuth
                and self.detector_azimuth == other.detector_azimuth
            )
        else:
            return False

    @property
    def reflection(self):
        return self.angleout > 0

    @property
    def anglenormin(self):
        """
        angle with surface normal (pointing inwards) in [0, 90]
        """
        return 90 - self.anglein

    @anglenormin.setter
    def anglenormin(self, value):
        self.anglein = 90 - value

    @property
    def anglenormout(self):
        """
        angle with surface normal (pointing inwards) in [0, 180]
        """
        return 90 + self.angleout

    @anglenormout.setter
    def anglenormout(self, value):
        self.angleout = value - 90

    @property
    def cosnormin(self):
        return np.cos(np.radians(self.anglenormin))

    @property
    def cosnormout(self):
        return np.cos(np.radians(self.anglenormout))

    @property
    def scatteringangle(self):
        return self.anglein + self.angleout

    @property
    def polar(self):
        return self.scatteringangle

    def xrayspectrumkwargs(self):
        return {
            "polar": np.radians(self.polar),
            "azimuth": np.radians(self.detector_azimuth),
        }

    def __str__(self):
        return "{}\nGeometry:\n In = {} deg\n Out = {} deg ({})\n Sample azimuth = {} deg\n Detector azimuth = {} deg".format(
            super(FlatSample, self).__str__(),
            self.anglein,
            self.angleout,
            "reflection" if self.reflection else "transmission",
            self.sample_azimuth,
            self.detector_azimuth,
        )

    def addtofisx(self, setup, cfg):
        # When self.angleout<0: works only for a single layer
        setup.setGeometry(self.anglein, abs(self.angleout))
        super(FlatSample, self).addtofisx(setup, cfg)


class SolidAngle(FlatSample):
    def __init__(self, solidangle=None, **kwargs):
        self.solidangle = solidangle  # srad
        super(SolidAngle, self).__init__(**kwargs)

    def __getstate__(self):
        state = super(SolidAngle, self).__getstate__()
        state["solidangle"] = self.solidangle
        return state

    def __setstate__(self, state):
        super(SolidAngle, self).__setstate__(state)
        self.solidangle = state["solidangle"]

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            if not (super(SolidAngle, self).__eq__(other)):
                return False
            return self.solidangle == other.solidangle
        else:
            return False

    def __str__(self):
        if self.solidangle is None:
            return super(SolidAngle, self).__str__()
        else:
            return "{}\n Solid angle = 4*pi*{} srad".format(
                super(SolidAngle, self).__str__(), self.solidangle / (4 * np.pi)
            )


class Centric(FlatSample):
    def __init__(self, distance=None, **kwargs):
        """
        Args:
            distance(num): distance (cm) to target
        """
        self.distance = distance
        super(Centric, self).__init__(**kwargs)

    def __getstate__(self):
        state = super(Centric, self).__getstate__()
        state["distance"] = self.distance
        return state

    def __setstate__(self, state):
        super(Centric, self).__setstate__(state)
        self.distance = state["distance"]

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            if not (super(Centric, self).__eq__(other)):
                return False
            return self.distance == other.distance
        else:
            return False

    @property
    def distance_rv(self):
        return self._distance

    @property
    def distance(self):
        return noisepropagation.E(self.distance_rv)

    @distance.setter
    def distance(self, value):
        if value is None:
            self._distance = None
        else:
            self._distance = units.Quantity(value, "cm")

    def calibrate_distance_manually(self, distance):
        self.distance = distance

    @property
    def solidangle_rv(self):
        return self.detector.solidangle_calc(
            activearea=self.detector.activearea_rv, distance=self.distance_rv
        )

    @property
    def solidangle(self):
        return noisepropagation.E(self.solidangle_rv)

    @solidangle.setter
    def solidangle(self, value):
        self.distance = self.detector.solidangle_calc(
            activearea=self.detector.activearea_rv, solidangle=value
        )

    def __str__(self):
        if self.distance is None:
            return super(Centric, self).__str__()
        else:
            if self.detector is None:
                solidangle = np.nan
            else:
                solidangle = "4*pi*{}".format(self.solidangle / (4 * np.pi))
            return "{}\n Distance = {:~}\n Solid angle = {} srad".format(
                super(Centric, self).__str__(), self.distance, solidangle
            )

    def addtopymca(self, setup, cfg):
        super(Centric, self).addtopymca(setup, cfg)
        cfg["concentrations"]["distance"] = self.distance.to("cm").magnitude

    def loadfrompymca(self, setup, cfg):
        super(Centric, self).loadfrompymca(setup, cfg)
        self.calibrate_distance_manually(cfg["concentrations"]["distance"])
