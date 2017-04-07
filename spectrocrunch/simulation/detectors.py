# -*- coding: utf-8 -*-
#
#   Copyright (C) 2017 European Synchrotron Radiation Facility, Grenoble, France
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

from six import with_metaclass

class DetectorMeta(type):
    """
    Metaclass used to register all detector classes inheriting from Detector
    """
    def __init__(cls, name, bases, dct):
        cls.registry[name.lower()] = cls
        super(DetectorMeta, cls).__init__(name, bases, dct)

class AreaDetector(with_metaclass(DetectorMeta, object)):
    """
    Class representing an area detector
    """
    registry = {}

    @classmethod
    def factory(cls, name):
        """
        Args:
            name(str): name of the detector

        Returns:
            AreaDetector
        """
        name = name.lower()
        if name in cls.registry:
            return cls.registry[name]()
        else:
            raise RuntimeError("Detector {} is not one of the registered detectors: {}".format(name, cls.registry.keys()))

    def __init__(self, etoadu=1, qe=1, aduoffset=0, darkcurrent=0, readoutnoise=0):
        """
        Args:
            etoadu(num): number of ADU per electron
            qe(num): detector quantum efficiency
            aduoffset(num): pixel intensity offset (ADU)
            darkcurrent(num): dark current (e/sec)
            readoutnoise(num): readout noise (e)
        """

        self.etoadu = float(etoadu)
        self.qe = float(qe)
        self.aduoffset = float(aduoffset)
        self.darkcurrent = float(darkcurrent)
        self.readoutnoise = float(readoutnoise)

class pcoedge55(AreaDetector):
    """
    PCO Edge 5.5
    """

    def __init__(self):
        super(pcoedge55, self).__init__(etoadu=1/0.45, qe=0.03, aduoffset=95.5, darkcurrent=7.4, readoutnoise=0.95)

