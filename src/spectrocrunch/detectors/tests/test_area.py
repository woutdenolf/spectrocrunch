import unittest
from .. import area
from ...patch import jsonpickle


class test_area(unittest.TestCase):
    def test_serialize(self):
        exclude = ("AreaDetector",)
        for name, cls in area.AreaDetector.clsregistry.items():
            if name not in exclude:
                d1 = cls()
                d2 = jsonpickle.loads(jsonpickle.dumps(d1))
                self.assertEqual(d1, d2)
