import unittest

from .. import visirlib
from ...patch.pint import ureg
from ...patch import jsonpickle

import numpy as np


class test_visirlib(unittest.TestCase):
    @unittest.skipIf(visirlib.PyTMM is None, "PyTMM not installed")
    def test_linatt(self):
        mat = visirlib.Material("main", "Al", "Rakic")
        linatt = mat.linear_attenuation_coefficient(
            ureg.Quantity([1.0332, 10.332], "micrometer")
        )
        linatt2 = np.asarray([9.8914e00, 8.8197e01])
        np.testing.assert_allclose(linatt, linatt2)

        mat = visirlib.Material("other", "air", "Ciddor")
        n = mat.refractive_index(ureg.Quantity([230.0, 960.0], "nm"))
        n2 = np.asarray([1.0003080029552, 1.0002742989071])
        np.testing.assert_allclose(n, n2)

    @unittest.skipIf(visirlib.PyTMM is None, "PyTMM not installed")
    def test_serialize(self):
        m1 = visirlib.Material("main", "Al", "Rakic")
        m2 = jsonpickle.loads(jsonpickle.dumps(m1))
        self.assertEqual(m1, m2)
