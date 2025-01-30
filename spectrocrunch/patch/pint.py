# -*- coding: utf-8 -*-

from __future__ import absolute_import
from pint import UnitRegistry

ureg = UnitRegistry()

ureg.define("classical_electron_radius = e^2/(4*pi*m_e*epsilon_0*c^2) = r_e")
ureg.define("percent = 1e-2*count = %")
ureg.define("permille = 1e-3*count = \u2030")
ureg.define("ppm = 1e-6*count")
ureg.define("ppb = 1e-9*count")
ureg.define("particles_per_mol = avogadro_number/mol")
