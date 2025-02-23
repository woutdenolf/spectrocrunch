import unittest
import numpy as np

from .. import emspectrum
from ...patch.pint import ureg
from ...patch import jsonpickle


class test_emspectrum(unittest.TestCase):
    def test_discrete(self):
        s1 = emspectrum.Discrete(ureg.Quantity([100, 200, 300], "nm"), [1, 1, 1])
        s2 = emspectrum.Discrete(ureg.Quantity([200, 400], "nm"), [1, 1])
        s3 = s1 + s2

        def func(x):
            return ureg.Quantity(x, "nm").to("keV", "spectroscopy")

        np.testing.assert_array_equal(s3.energies, func([100, 200, 300, 400]))
        np.testing.assert_array_equal(s3.ratios, [0.2, 0.4, 0.2, 0.2])
        s4 = emspectrum.Discrete(ureg.Quantity([150, 350], "nm"), [1, 2])
        self.assertEqual(s4.sample(s3), 1.5 * 1 / 3.0 + 1 * 2 / 3.0)

    def test_dark(self):
        s1 = emspectrum.Dark()
        s4 = emspectrum.Discrete(ureg.Quantity([150, 350], "nm"), [1, 2])
        self.assertEqual(s4.sample(s1), 0)

    def test_serialize(self):
        s1 = emspectrum.Discrete(ureg.Quantity([100, 300], "nm"), [1, 1])
        s2 = jsonpickle.loads(jsonpickle.dumps(s1))
        self.assertEqual(s1, s2)
        s1 = emspectrum.Dark()
        s2 = jsonpickle.loads(jsonpickle.dumps(s1))
        self.assertEqual(s1, s2)
