# -*- coding: utf-8 -*-

import unittest

from .. import ft

import numpy as np


class test_ft(unittest.TestCase):
    """http://dsp.stackexchange.com/questions/633/what-data-should-i-use-to-test-an-fft-implementation-and-what-accuracy-should-i
       http://citeseer.ist.psu.edu/viewdoc/summary?doi=10.1.1.38.3924
    """

    def test_dirac_delta(self):
        for convention in [ft.fftConvention.numpy, ft.fftConvention.idl]:
            for dx in [1, 0.5]:
                for dy in [1, 0.5]:
                    for nx in range(4, 6):
                        sig = np.zeros(nx)
                        sig[0] = 1
                        ftsig1 = ft.fft(sig, dx=dx, normconvention=convention)
                        ftsig2 = np.ones(nx, dtype=np.complex)
                        if convention == ft.fftConvention.idl:
                            ftsig2 /= nx
                        np.testing.assert_allclose(ftsig1, ftsig2)

                        ftsig = np.ones(nx)
                        sig1 = ft.ifft(ftsig, dx=dx, normconvention=convention)
                        sig2 = np.zeros(nx, dtype=np.complex)
                        sig2[0] = 1
                        if convention == ft.fftConvention.idl:
                            sig2[0] = nx
                        np.testing.assert_allclose(sig1, sig2, atol=1e-10)

                        for ny in range(4, 6):
                            sig = np.zeros((ny, nx))
                            sig[0, 0] = 1
                            ftsig1 = ft.fft2(
                                sig, dx=dx, dy=dy, normconvention=convention
                            )
                            ftsig2 = np.ones((ny, nx), dtype=np.complex)
                            if convention == ft.fftConvention.idl:
                                ftsig2 /= nx * ny
                            np.testing.assert_allclose(ftsig1, ftsig2)

                            ftsig = np.ones((ny, nx))
                            sig1 = ft.ifft2(
                                ftsig, dx=dx, dy=dy, normconvention=convention
                            )
                            sig2 = np.zeros((ny, nx), dtype=np.complex)
                            sig2[0, 0] = 1
                            if convention == ft.fftConvention.idl:
                                sig2[0, 0] = nx * ny
                            np.testing.assert_allclose(sig1, sig2, atol=1e-10)

    def test_consistent(self):
        for convention in [ft.fftConvention.numpy, ft.fftConvention.idl]:
            for dx in [1, 0.5]:
                for dy in [1, 0.5]:
                    for nx in range(4, 6):
                        sig = np.random.rand(nx)
                        np.testing.assert_allclose(
                            sig,
                            ft.ifft(
                                ft.fft(sig, dx=dx, normconvention=convention),
                                dx=dx,
                                normconvention=convention,
                            ),
                        )
                        np.testing.assert_allclose(
                            sig,
                            ft.fft(
                                ft.ifft(sig, dx=dx, normconvention=convention),
                                dx=dx,
                                normconvention=convention,
                            ),
                        )
                        np.testing.assert_allclose(
                            ft.fft(sig, dx=dx), ft.fft(sig, dx=dx, u=ft.fftfreq(nx, dx))
                        )
                        for ny in range(4, 6):
                            sig = np.random.rand(ny, nx)
                            np.testing.assert_allclose(
                                sig,
                                ft.ifft2(
                                    ft.fft2(
                                        sig, dx=dx, dy=dy, normconvention=convention
                                    ),
                                    dx=dx,
                                    dy=dy,
                                    normconvention=convention,
                                ),
                            )
                            np.testing.assert_allclose(
                                sig,
                                ft.fft2(
                                    ft.ifft2(
                                        sig, dx=dx, dy=dy, normconvention=convention
                                    ),
                                    dx=dx,
                                    dy=dy,
                                    normconvention=convention,
                                ),
                            )
                            np.testing.assert_allclose(
                                ft.fft2(sig, dx=dx, dy=dy),
                                ft.fft2(
                                    sig,
                                    dx=dx,
                                    dy=dy,
                                    u=ft.fftfreq(nx, dx),
                                    v=ft.fftfreq(ny, dy),
                                ),
                            )

    def test_linear(self):
        for dx in [1, 0.5]:
            for dy in [1, 0.5]:
                for nx in range(4, 6):
                    x = np.random.rand(nx)
                    y = np.random.rand(nx)
                    a = np.random.rand(1)[0]
                    b = np.random.rand(1)[0]
                    np.testing.assert_allclose(
                        ft.fft(a * x + b * y, dx=dx),
                        a * ft.fft(x, dx=dx) + b * ft.fft(y, dx=dx),
                    )
                    np.testing.assert_allclose(
                        ft.ifft(a * x + b * y, dx=dx),
                        a * ft.ifft(x, dx=dx) + b * ft.ifft(y, dx=dx),
                    )
                    for ny in range(4, 6):
                        x = np.random.rand(ny, nx)
                        y = np.random.rand(ny, nx)
                        a = np.random.rand(1)[0]
                        b = np.random.rand(1)[0]
                        np.testing.assert_allclose(
                            ft.fft2(a * x + b * y, dx=dx),
                            a * ft.fft2(x, dx=dx) + b * ft.fft2(y, dx=dx),
                        )
                        np.testing.assert_allclose(
                            ft.ifft2(a * x + b * y, dx=dx),
                            a * ft.ifft2(x, dx=dx) + b * ft.ifft2(y, dx=dx),
                        )

    def test_comp_np(self):
        for nx in range(4, 6):
            sig = np.random.rand(nx)
            np.testing.assert_allclose(
                np.fft.fft(sig), ft.fft(sig, normconvention=ft.fftConvention.numpy)
            )
            np.testing.assert_allclose(
                np.fft.ifft(sig), ft.ifft(sig, normconvention=ft.fftConvention.numpy)
            )
            for ny in range(4, 6):
                sig = np.random.rand(ny, nx)
                np.testing.assert_allclose(
                    np.fft.fft2(sig),
                    ft.fft2(sig, normconvention=ft.fftConvention.numpy),
                )
                np.testing.assert_allclose(
                    np.fft.ifft2(sig),
                    ft.ifft2(sig, normconvention=ft.fftConvention.numpy),
                )

        # Force using non-numpy code
        for nx in range(4, 6):
            sig = np.random.rand(nx)
            u = ft.fftfreq(nx)
            np.testing.assert_allclose(
                np.fft.fft(sig), ft.fft(sig, normconvention=ft.fftConvention.numpy, u=u)
            )
            np.testing.assert_allclose(
                np.fft.ifft(sig),
                ft.ifft(sig, normconvention=ft.fftConvention.numpy, u=u),
            )
            for ny in range(4, 6):
                v = ft.fftfreq(ny)
                sig = np.random.rand(ny, nx)
                np.testing.assert_allclose(
                    np.fft.fft2(sig),
                    ft.fft2(sig, normconvention=ft.fftConvention.numpy, u=u, v=v),
                )
                np.testing.assert_allclose(
                    np.fft.ifft2(sig),
                    ft.ifft2(sig, normconvention=ft.fftConvention.numpy, u=u, v=v),
                )

    def test_shift(self):
        for dx in [1, 0.5]:
            for dy in [1, 0.5]:
                for nx in range(4, 6):
                    sig1 = np.zeros(nx)
                    sig1[0] = 1
                    sig2 = np.zeros(nx)
                    ox = 1 + np.random.rand(1)[0] * (nx - 2)
                    ox = int(ox)
                    sig2[ox] = 1
                    ftsig1 = ft.fft(sig1, dx=dx)
                    ftsig2 = ft.fft(sig2, dx=dx) * np.exp(
                        2j * np.pi * ft.fftfreq(nx, d=dx) * ox * dx
                    )
                    np.testing.assert_allclose(ftsig1, ftsig2)

                    for ny in range(4, 6):
                        sig1 = np.zeros((ny, nx))
                        sig1[0, 0] = 1
                        sig2 = np.zeros((ny, nx))
                        ox, oy = 1 + np.random.rand(2) * (np.array([nx, ny]) - 2)
                        ox = int(ox)
                        oy = int(oy)
                        sig2[oy, ox] = 1
                        ftsig1 = ft.fft2(sig1, dx=dx, dy=dy)
                        ftsig2 = ft.fft2(sig2, dx=dx, dy=dy) * np.exp(
                            2j
                            * np.pi
                            * (
                                np.add.outer(
                                    ft.fftfreq(ny, d=dy) * oy * dy,
                                    ft.fftfreq(nx, d=dx) * ox * dx,
                                )
                            )
                        )
                        np.testing.assert_allclose(ftsig1, ftsig2)

    def test_freq(self):
        for freqconvention in [ft.fftConvention.numpy, ft.fftConvention.idl]:
            for dx in [1, 0.5]:
                for dy in [1, 0.5]:
                    for nx in range(4, 6):
                        sig = np.random.rand(nx)
                        ftsig = ft.fftshift(
                            ft.fft(sig, dx=dx), freqconvention=freqconvention
                        )
                        freq = ft.fftshift(
                            ft.fftfreq(nx, d=dx, freqconvention=freqconvention),
                            freqconvention=freqconvention,
                        )
                        self.assertTrue(all(np.diff(freq) > 0))
                        if nx % 2 == 1:
                            self.assertEqual(ftsig[0], np.conj(ftsig[-1]))

                        for ny in range(4, 6):
                            sig = np.random.rand(ny, nx)
                            ftsig = ft.fftshift(
                                ft.fft2(sig, dx=dx, dy=dy),
                                freqconvention=freqconvention,
                            )
                            u = ft.fftshift(
                                ft.fftfreq(nx, d=dx, freqconvention=freqconvention),
                                freqconvention=freqconvention,
                            )
                            v = ft.fftshift(
                                ft.fftfreq(ny, d=dy, freqconvention=freqconvention),
                                freqconvention=freqconvention,
                            )
                            self.assertTrue(all(np.diff(u) > 0))
                            self.assertTrue(all(np.diff(v) > 0))
                            if nx % 2 == 1 and ny % 2 == 1:
                                self.assertAlmostEqual(
                                    ftsig[0, 0], np.conj(ftsig[-1, -1])
                                )

    def test_subregion(self):
        for freqconvention in [ft.fftConvention.numpy, ft.fftConvention.idl]:
            for dx in [1, 0.5]:
                for dy in [1, 0.5]:
                    for nx in range(4, 6):
                        sig = np.random.rand(nx)
                        ftsig = ft.fft(sig, dx=dx)
                        # subregion in real space:
                        np.testing.assert_allclose(
                            sig[1 : nx - 1], ft.ifft(ftsig, dx=dx, x0=1, x1=nx - 2)
                        )
                        # subregion in Fourier space:
                        u = ft.fftshift(
                            ft.fftfreq(nx, d=dx, freqconvention=freqconvention),
                            freqconvention=freqconvention,
                        )
                        np.testing.assert_allclose(
                            ft.fftshift(ftsig, freqconvention=freqconvention)[
                                1 : nx - 1
                            ],
                            ft.fft(sig, dx=dx, u=u[1 : nx - 1]),
                        )
                        for ny in range(4, 6):
                            sig = np.random.rand(ny, nx)
                            ftsig = ft.fft2(sig, dx=dx, dy=dy)
                            # subregion in real space:
                            np.testing.assert_allclose(
                                sig[1 : ny - 1, 1 : nx - 1],
                                ft.ifft2(
                                    ftsig,
                                    dx=dx,
                                    x0=1,
                                    x1=nx - 2,
                                    dy=dy,
                                    y0=1,
                                    y1=ny - 2,
                                ),
                            )
                            # subregion in Fourier space:
                            u = ft.fftshift(
                                ft.fftfreq(nx, d=dx, freqconvention=freqconvention),
                                freqconvention=freqconvention,
                            )
                            v = ft.fftshift(
                                ft.fftfreq(ny, d=dy, freqconvention=freqconvention),
                                freqconvention=freqconvention,
                            )
                            np.testing.assert_allclose(
                                ft.fftshift(ftsig, freqconvention=freqconvention)[
                                    1 : ny - 1, 1 : nx - 1
                                ],
                                ft.fft2(
                                    sig, dx=dx, u=u[1 : nx - 1], dy=dy, v=v[1 : ny - 1]
                                ),
                            )

    def test_center(self):
        for freqconvention in [ft.fftConvention.numpy, ft.fftConvention.idl]:
            for dx in [1, 0.5]:
                for dy in [1, 0.5]:
                    for nx in range(4, 6):
                        sig = np.random.rand(nx)
                        np.testing.assert_allclose(
                            ft.fftshift(ft.fft(sig)), ft.fft(sig, centered=True)
                        )
                        np.testing.assert_allclose(
                            ft.fftshift(ft.fftfreq(nx)), ft.fftfreq(nx, centered=True)
                        )
                        np.testing.assert_allclose(
                            sig, ft.ifft(ft.fft(sig, centered=True), centered=True)
                        )
                        for ny in range(4, 6):
                            sig = np.random.rand(ny, nx)
                            np.testing.assert_allclose(
                                ft.fftshift(ft.fft2(sig)), ft.fft2(sig, centered=True)
                            )
                            np.testing.assert_allclose(
                                sig,
                                ft.ifft2(ft.fft2(sig, centered=True), centered=True),
                            )


def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_ft("test_comp_np"))
    testSuite.addTest(test_ft("test_consistent"))
    testSuite.addTest(test_ft("test_dirac_delta"))
    testSuite.addTest(test_ft("test_linear"))
    testSuite.addTest(test_ft("test_shift"))
    testSuite.addTest(test_ft("test_freq"))
    testSuite.addTest(test_ft("test_subregion"))
    testSuite.addTest(test_ft("test_center"))
    return testSuite


if __name__ == "__main__":
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
