# -*- coding: utf-8 -*-

import unittest

from .. import spec


class test_spec(unittest.TestCase):

    def test_cmd_parser(self):
        cmds = [("zapimage sampy 0 1 10 sampz 2 3 11 100",
                 {"time": lambda x: x.to("ms").magnitude == 100}),
                ("zapimage sampy 0 1 10 100 sampz 2 3 11",
                 {"time": lambda x: x.to("ms").magnitude == 100}),
                ("puzzle sampy 0 1 10 sampz 2 3 11 100",
                 {"time": lambda x: x.to("ms").magnitude == 100}),
                ("mesh sampy 0 1 10 sampz 2 3 11 0.1",
                 {"time": lambda x: x.to("ms").magnitude == 100}),
                ("zapline sampy 0 1 10 100",
                 {"time": lambda x: x.to("ms").magnitude == 100}),
                ("ascan sampy 0 1 10 0.1",
                 {"time": lambda x: x.to("ms").magnitude == 100}),
                ("zapenergy SUM 10 100",
                 {"time": lambda x: x.to("ms").magnitude == 100}),
                ("zapenergy SUM2 10 100",
                 {"time": lambda x: x.to("ms").magnitude == 100}),
                ("invalid", {"name": lambda x: x == "unknown"})]

        p = spec.cmd_parser()
        for cmd, checks in cmds:
            r = p.parse(cmd)
            for k, func in checks.items():
                self.assertTrue(func(r[k]))


def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_spec("test_cmd_parser"))
    return testSuite


if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
