import unittest


def main_test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    return testSuite


if __name__ == "__main__":
    import sys

    mysuite = main_test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
