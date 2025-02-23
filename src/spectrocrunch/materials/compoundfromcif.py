import os

from . import compound
from . import types
from ..resources import resource_filename

try:
    import iotbx.cif as iotbxcif
except ImportError:
    iotbxcif = None


class CompoundFromCif(compound.Compound):
    """Interface to a compound defined by a cif file"""

    def __init__(self, filename, density=0, name=None):
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
            raise IOError("Cif file %s not found." % filename)
        self.ciffile = f

        # cctbx.xray.structure.structure
        self.structure = (
            iotbxcif.reader(file_path=f).build_crystal_structures().values()[0]
        )

        # Extract necessary information
        scatterers = {}
        for s in self.structure.scatterers():
            e = s.scattering_type
            if e in scatterers:
                scatterers[e] += s.occupancy * s.multiplicity()
            else:
                scatterers[e] = s.occupancy * s.multiplicity()

        if density == 0:
            density = self.structure.crystal_density()

        super(CompoundFromCif, self).__init__(
            scatterers.keys(),
            scatterers.values(),
            types.fraction.mole,
            density,
            name=name,
        )

    def _get_cif_name(self, filename):
        """Get file from the database if it doesn't exist"""
        name = resource_filename(filename)
        if os.path.isfile(name):
            return name
        else:
            f = os.path.join(
                os.path.dirname(os.path.abspath(__file__)),
                "cif",
                os.path.splitext(os.path.basename(name))[0] + ".cif",
            )
            if os.path.isfile(f):
                return f
            else:
                return None
