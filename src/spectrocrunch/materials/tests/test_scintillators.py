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


def main_test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_scintillators("test_serialize"))
    return testSuite


if __name__ == "__main__":
    import sys

    mysuite = main_test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
