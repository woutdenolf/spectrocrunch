import unittest

from .. import stoichiometry

import numpy as np


class test_stoichiometry(unittest.TestCase):
    def test_fraction(self):
        n = np.asarray([0.1, 0.6, 0.3])
        MM = np.asarray([5.6, 3.5, 7.9])
        rho = 5.3

        x = stoichiometry.frac_mole_to_weight(n, MM)
        n2 = stoichiometry.frac_weight_to_mole(x, MM)
        np.testing.assert_allclose(n, n2)

        x = stoichiometry.frac_mole_to_volume(n, rho, MM)
        n2 = stoichiometry.frac_volume_to_mole(x, rho, MM)
        np.testing.assert_allclose(n, n2)

        x = stoichiometry.frac_weight_to_volume(n, rho)
        n2 = stoichiometry.frac_volume_to_weight(x, rho)
        np.testing.assert_allclose(n, n2)
