import unittest
import os

from ...io import nxfs
from ...patch import xraylib
from .xrfmap import XrfMapGenerator
from .test_task import test_task


class test_task_xrf(test_task):

    @unittest.skipIf(xraylib.XRayInit is None, "xraylib not installed")
    def test_nxqxrf(self):
        fixed = {'plot': False,
                 'dark': False,
                 'time': 1,
                 'gaindiodeI0': 1e8,
                 'gaindiodeIt': 1e7
        }
        variable = [{'I0_counts': 300, 'It_counts': 30,
                     'dark': True},
                    {'I0_counts': 400000, 'It_counts': 100000,
                     'energy': 7},
                    {'I0_counts': 200000, 'It_counts': 100000,
                     'energy': 7.1}
        ]
        h5filename = os.path.join(self.dir.path, 'test.h5')
        outputparent = nxfs.Path('/', h5file=h5filename)['entry']
        parameters = {
            'method': 'xrfgeometry',
            'outputparent': outputparent,
            'geometry': 'sxm',
            'init': {},
            'fixed': fixed,
            'variable': variable
        }
        proc2 = self._run_task(parameters, None)
        parameters['geometry'] = 'SXM'
        proc3 = self._run_task(parameters, None)
        self._check_reproc(proc2, proc3)

    @unittest.skipIf(xraylib.XRayInit is None, "xraylib not installed")
    @unittest.skip('not implemented yet')
    def test_nxxiaedf(self):
        xrfmap = XrfMapGenerator(nmaps=1)
        xrfmap.generate(self.dir.path, 'test')

        # Needs a quantification run
        h5filename = os.path.join(self.dir.path, 'geometries.h5')
        seldetectors = [(det,) for det in range(xrfmap.ndet)]
        parameters = xrfmap.taskparams_geometry(seldetectors)
        parameters['method'] = 'xrfgeometry'
        parameters['outputparent'] = h5filename + '::/xrf'
        proc1 = self._run_task(parameters, None)

        # Output is an NXentry (not an NXprocess)
        h5filename = os.path.join(self.dir.path, 'test.h5')
        parameters = {
            'method': 'xiaedftonx',
            'outputparent': h5filename,
            'path': xrfmap.path,
            'radix': xrfmap.radix,
            'number': xrfmap.scannumbers[0],
            'instrument': xrfmap.instrument.lower()
        }
        parameters.update(xrfmap.taskparams_pymca(seldetectors))

        proc2 = self._run_task(parameters, proc1, outputnxprocess=False)
        parameters['instrument'] = xrfmap.instrument.upper()
        proc3 = self._run_task(parameters, proc1, outputnxprocess=False)
        self._check_reproc(proc2, proc3)


def test_suite():
    '''Test suite including all test suites'''
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_task_xrf('test_nxqxrf'))
    testSuite.addTest(test_task_xrf('test_nxxiaedf'))
    return testSuite


if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
