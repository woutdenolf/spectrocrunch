import unittest

from .. import randomdata


class test_randomdata(unittest.TestCase):
    def test_native(self):
        for _ in range(500):
            o = randomdata.factory(types=("native",))
            # Check equality (shuffles unsorted types)
            self.assertEqual(o, o)

    def test_all(self):
        for _ in range(500):
            o = randomdata.factory()
            self.assertTrue("random" in str(o).lower())
            # Make sure raw data is generated
            data = o.data
            try:
                us = str(data)
            except UnicodeDecodeError:
                us = data.decode("latin1")
            self.assertFalse("random" in us.lower())
            # Check equality (shuffles unsorted types)
            self.assertEqual(o, o)
