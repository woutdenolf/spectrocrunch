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


def main_test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_area("test_serialize"))
    return testSuite


if __name__ == "__main__":
    import sys

    mysuite = main_test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
