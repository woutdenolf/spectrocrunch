import unittest
from .. import listtools


class test_listtools(unittest.TestCase):
    def test_unique(self):
        a, b = [2, 0, 0, 1, 0, 1], [0, 1, 2, 3, 4, 5]
        x, y = listtools.unique2lists(a, b, add=False)
        self.assertEqual(x, [2, 0, 1])
        self.assertEqual(y, [0, 1, 3])
        x, y = listtools.unique2lists(a, b, add=True)
        self.assertEqual(x, [2, 0, 1])
        self.assertEqual(y, [0, 7, 8])

    def test_sort(self):
        a, b = [2, 0, 0, 1, 0, 1], [0, 1, 2, 3, 4, 5]
        x, y = listtools.sort2lists(a, b)
        self.assertEqual(x, [0, 0, 0, 1, 1, 2])
        self.assertEqual(y, [1, 2, 4, 3, 5, 0])


def main_test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_listtools("test_unique"))
    testSuite.addTest(test_listtools("test_sort"))
    return testSuite


if __name__ == "__main__":
    import sys

    mysuite = main_test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
