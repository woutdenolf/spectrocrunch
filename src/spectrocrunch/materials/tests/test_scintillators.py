import unittest

from .. import scintillators
from ...patch import jsonpickle


class test_scintillators(unittest.TestCase):
    @unittest.skipIf(
        scintillators.compound.element.xraylib.XRayInit is None, "xraylib not installed"
    )
    def test_serialize(self):
        exclude = ()
        for name, cls in scintillators.Scintillator.clsregistry.items():
            if name not in exclude:
                s1 = cls()
                s2 = jsonpickle.loads(jsonpickle.dumps(s1))
                self.assertEqual(s1, s2)
