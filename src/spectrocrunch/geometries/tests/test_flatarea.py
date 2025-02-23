import unittest

from .. import flatarea
from ...patch import jsonpickle


class test_flatarea(unittest.TestCase):
    def test_serialize(self):
        exclude = ()
        for name, cls in flatarea.Geometry.clsregistry.items():
            if name not in exclude:
                g1 = cls()
                g2 = jsonpickle.loads(jsonpickle.dumps(g1))
                self.assertEqual(g1, g2)
