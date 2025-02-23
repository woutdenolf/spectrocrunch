import numpy as np

from .. import xray
from ...utils import units
from ...patch import jsonpickle
from ...materials import compoundfromname
from ...testutils.subtest import TestCase


class test_xray(TestCase):
    def test_interpolate(self):
        o1 = xray.KB(kind="linear")
        o1.set_transmission(7, 0.2)
        o1.set_transmission(units.Quantity(7400, "eV"), 0.8)

        def transmission(x):
            return o1.transmission(units.Quantity(x, "keV"))

        self.assertEqual(transmission(7.2), 0.5)
        np.testing.assert_allclose(transmission([7.2]), [0.5])
        np.testing.assert_allclose(transmission([7.1, 7.2, 7.3]), [0.35, 0.5, 0.65])

    def test_serialize(self):
        with self.skipContext():
            exclude = ()
            for name, cls in xray.XrayOptics.clsregistry.items():
                if name not in exclude:
                    with self.subTest(name=name):
                        if name == "Filter":
                            if compoundfromname.xraylib is None:
                                self.skipTest("xraylib not installed")
                            material = compoundfromname.compoundfromname("air")
                            o1 = cls(material=material, thickness=1)
                        else:
                            o1 = cls()
                        o2 = jsonpickle.loads(jsonpickle.dumps(o1))
                        self.assertEqual(o1, o2)
