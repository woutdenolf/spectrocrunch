# -*- coding: utf-8 -*-

import unittest

from ..types import transformationType
from . import helper_teststack

import numpy as np
import matplotlib.pyplot as plt


class test_teststack(unittest.TestCase):
    def test_data(self):
        types = [
            transformationType.translation,
            transformationType.rigid,
            transformationType.similarity,
            transformationType.affine,
            transformationType.projective,
        ]
        for t in types:
            if t == transformationType.translation:
                lst = [True, False]
            else:
                lst = [False]

            for vector in lst:
                for transposed in lst:
                    listofstacks, COFrelative, stackdim = helper_teststack.data(
                        t, vector=vector, transposed=transposed
                    )
                    self.assertIsInstance(listofstacks, list)
                    self.assertIsInstance(listofstacks[0], np.ndarray)
                    self.assertEqual(len(listofstacks[0].shape), 3)
                    self.assertTrue(
                        all(s.shape == listofstacks[0].shape for s in listofstacks)
                    )

                    if False:
                        for s in listofstacks:
                            for i in range(s.shape[stackdim]):
                                if stackdim == 0:
                                    self.plot(s[i, ...])
                                elif stackdim == 1:
                                    self.plot(s[:, i, :])
                                else:
                                    self.plot(s[..., i])
                            break  # show only one

    def plot(self, img):
        plt.figure(1)
        plt.subplot(111)
        plt.imshow(img, origin="lower", interpolation="nearest")
        plt.pause(0.1)


def main_test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_teststack("test_data"))
    return testSuite


if __name__ == "__main__":
    import sys

    mysuite = main_test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
