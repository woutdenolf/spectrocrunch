# -*- coding: utf-8 -*-

execfile("initcctbx.py")

# Initialize CCTBX environment as . initcctbx

import iotbx.cif

if __name__ == "__main__":

    ciffile = "/data/id21/inhouse/wout/dev/SpectroCrunch/spectrocrunch/resources/cif/calcite.cif"
    struct = iotbx.cif.reader(file_path=ciffile).build_crystal_structures().values()[0]
    scat = struct.scatterers()
    for s in scat:
        print(s.occupancy * s.multiplicity())
        print(s.scattering_type)
