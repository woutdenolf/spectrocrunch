import unittest

from .. import base
from ...patch import jsonpickle
from ...materials import compoundfromname
from ...materials import element
from ...materials import mixture


class test_base(unittest.TestCase):
    @unittest.skipIf(compoundfromname.xraylib is None, "xraylib not installed")
    def test_serialize(self):
        ca = element.Element("Ca")
        calcite = compoundfromname.compoundfromname("calcite")
        hematite = compoundfromname.compoundfromname("hematite")
        goethite = compoundfromname.compoundfromname("goethite")
        mix = mixture.Mixture(
            [ca, hematite, goethite, calcite],
            [0.25, 0.25, 0.25, 0.25],
            name="iron oxides",
        )
        attenuators = {"A": ca, "B": calcite, "C": mix}

        d1 = base.Material()
        d2 = jsonpickle.loads(jsonpickle.dumps(d1))
        self.assertEqual(d1, d2)
        d1 = base.Material(attenuators=attenuators)
        d2 = jsonpickle.loads(jsonpickle.dumps(d1))
        self.assertEqual(d1, d2)

        d1 = base.SolidState()
        d2 = jsonpickle.loads(jsonpickle.dumps(d1))
        self.assertEqual(d1, d2)
        d1 = base.SolidState(attenuators=attenuators, ehole=3.6)
        d2 = jsonpickle.loads(jsonpickle.dumps(d1))
        self.assertEqual(d1, d2)

        d1 = base.CentricCone()
        d2 = jsonpickle.loads(jsonpickle.dumps(d1))
        self.assertEqual(d1, d2)
        d1 = base.CentricCone(attenuators=attenuators, ehole=3.6, activearea=0.8)
        d2 = jsonpickle.loads(jsonpickle.dumps(d1))
        self.assertEqual(d1, d2)
