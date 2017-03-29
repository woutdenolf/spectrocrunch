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

from .compound import compound
from .types import fractionType
try:
    import iotbx.cif as iotbxcif
except ImportError:
    iotbxcif = None
import os

class compoundfromcif(compound):
    """Interface to a compound defined by a cif file
    """

    def __init__(self,filename,density=0,name=None):
        """
        Args:
            filename(str): cif file name
            density(Optional[float]): compound density
            name(Optional[str]): compound name
        Raises:
            IOError: If the file doesn't exist
            RuntimeError: 
        """
        if iotbxcif is None:
            raise RuntimeError("cctbx is required to read cif files")

        f = self._get_cif_name(filename)
        if f is None:
            raise IOError("Cif file %s not found."%filename)
        self.ciffile = f

        # cctbx.xray.structure.structure
        self.structure = iotbxcif.reader(file_path=f).build_crystal_structures().values()[0]

        # Extract necessary information
        scatterers = {}
        for s in self.structure.scatterers():
            e = s.scattering_type
            if e in scatterers:
                scatterers[e] += s.occupancy*s.multiplicity()
            else:
                scatterers[e] = s.occupancy*s.multiplicity()

        if density==0:
            density = self.structure.crystal_density()

        super(compoundfromcif,self).__init__(scatterers.keys(),scatterers.values(),fractionType.mole,density,name=name)

    def _get_cif_name(self,name):
        """Get file from the database if it doesn't exist
        """
        if os.path.isfile(name):
            return name
        else:
            f = os.path.join(os.path.dirname(os.path.abspath(__file__)),"cif",os.path.splitext(os.path.basename(name))[0]+".cif")
            if os.path.isfile(f):
                return f
            else:
                return None
