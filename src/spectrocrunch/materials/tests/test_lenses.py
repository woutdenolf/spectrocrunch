import unittest

from .. import lenses
from ...patch import jsonpickle


class test_lenses(unittest.TestCase):
    @unittest.skipIf(lenses.visirlib.PyTMM is None, "PyTMM not installed")
    def test_serialize(self):
        exclude = ()
        for name, cls in lenses.Lens.clsregistry.items():
            if name not in exclude:
                l1 = cls()
                l2 = jsonpickle.loads(jsonpickle.dumps(l1))
                self.assertEqual(l1, l2)
