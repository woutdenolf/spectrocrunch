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

import os
import warnings
from glob import glob
try:
    import PyTMM
    import PyTMM.refractiveIndex
except ImportError:
    PyTMM = None
    warnings.warn("PyTMM is not installed", ImportWarning)
else:
    from ..utils.instance import isarray
    from ..patch.pint import ureg

    root = os.path.dirname(os.path.realpath(PyTMM.__file__))
    files = glob(os.path.join(root, '*', 'library.yml'))
    if files:
        path = os.path.dirname(files[0])
        db = PyTMM.refractiveIndex.RefractiveIndex(path)
    else:
        db = PyTMM.refractiveIndex.RefractiveIndex()

    class Material(PyTMM.refractiveIndex.Material):

        def __init__(self, shelf, book, page):
            shelf = shelf.lower()
            self._material_id = shelf, book, page
            filename = db.getMaterialFilename(shelf, book, page)
            PyTMM.refractiveIndex.Material.__init__(self, filename)

        def __getstate__(self):
            return {'_material_id': self._material_id}

        def __setstate__(self, state):
            material_id = state['_material_id']
            o = self.__class__(*material_id)
            self._material_id = material_id
            self.refractiveIndex = o.refractiveIndex
            self.extinctionCoefficient = o.extinctionCoefficient

        def __eq__(self, other):
            if isinstance(other, self.__class__):
                return self._material_id == other._material_id
            else:
                return False

        def __ne__(self, other):
            return not self.__eq__(other)

        def linear_attenuation_coefficient(self, lines):
            """Linear absorption coefficient (1/cm)

            Args:
                lines(ureg.Quantity): keV, nm, ...

            Returns:
                num|array
            """
            wl = lines.to("nm", "spectroscopy").magnitude
            if isarray(wl):
                return [self.getExtinctionCoefficient(l) for l in wl]
            else:
                return self.getExtinctionCoefficient(wl)

        def refractive_index(self, lines):
            """Refractive index

            Args:
                lines(array(lines)): keV, nm, ...

            Returns:
                num|array
            """
            wl = lines.to("nm", "spectroscopy").magnitude
            if isarray(wl):
                return [self.getRefractiveIndex(l) for l in wl]
            else:
                return self.getRefractiveIndex(wl)
