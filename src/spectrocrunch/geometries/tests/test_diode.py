import unittest

from .. import diode
from ...patch import jsonpickle


class test_diode(unittest.TestCase):
    def test_serialize(self):
        exclude = ("DiodeGeometry",)
        for name, cls in diode.DiodeGeometry.clsregistry.items():
            if name not in exclude:
                d1 = cls()
                d2 = jsonpickle.loads(jsonpickle.dumps(d1))
                self.assertEqual(d1, d2)
