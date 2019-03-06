# -*- coding: utf-8 -*-
#
#   Copyright (C) 2018 European Synchrotron Radiation Facility, Grenoble, France
#
#   Principal author:   Wout De Nolf (wout.de_nolf@esrf.eu)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

# -*- coding: utf-8 -*-
#
#   Copyright (C) 2018 European Synchrotron Radiation Facility, Grenoble, France
#
#   Principal author:   Wout De Nolf (wout.de_nolf@esrf.eu)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import unittest
import numpy as np
import itertools

from .. import interpolate


class test_interpolate(unittest.TestCase):

    def _generate(self, shape):
        axes = [range(n) for n in shape]
        data = np.random.normal(size=shape).astype(np.float32)*100
        return data, axes

    def test_nd(self):
        def identity(x): return x
        def meshgrid(axes): return np.meshgrid(*axes, indexing='ij')

        parameters = ((0, 1, 2), (True,),
                      ((interpolate.interpolate_regular, identity),
                       (interpolate.interpolate_irregular, meshgrid)))

        for parameters in itertools.product(*parameters):
            degree, asgrid, (interp, meshgrid) = parameters

            rtol = 1e-3

            # 1D
            data1, axes = self._generate((6,))
            axnew = (0,)
            data2 = interp(data1, axes, axnew, degree=degree, asgrid=asgrid)
            data3 = data1[axnew]
            np.testing.assert_allclose(data2, data3, rtol=rtol)
            axnew = ([0, 5],)
            data2 = interp(data1, axes, axnew, degree=degree, asgrid=asgrid)
            data3 = data1[axnew]
            np.testing.assert_allclose(data2, data3, rtol=rtol)

            # 2D
            data1, axes = self._generate((6, 8))
            if asgrid:
                axnew = ([0, 3, 5], [1, 3, 6])
            else:
                axnew = ([0, 3, 5], [1, 3])
            data2 = interp(data1, meshgrid(axes), axnew,
                           degree=degree, asgrid=asgrid)
            if asgrid:
                axnew = tuple(np.meshgrid(*axnew, indexing='ij'))
            data3 = data1[axnew]
            np.testing.assert_allclose(data2, data3, rtol=rtol)

            # 3D
            data1, axes = self._generate((6, 8, 9))
            if asgrid:
                axnew = ([0, 3, 5], [1, 3, 6], [2, 5, 7])
            else:
                axnew = ([0, 3, 5], [1, 3], [2, 7])
            data2 = interp(data1, meshgrid(axes), axnew,
                           degree=degree, asgrid=asgrid)
            if asgrid:
                axnew = tuple(np.meshgrid(*axnew, indexing='ij'))
            data3 = data1[axnew]
            np.testing.assert_allclose(data2, data3, rtol=rtol)


def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_interpolate("test_nd"))
    return testSuite


if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
