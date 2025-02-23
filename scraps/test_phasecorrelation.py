import numpy as np
import spectrocrunch.math.ft as ft

centered = False

n = 5
sig1 = np.zeros(n)
sig2 = np.zeros(n)

m = 3
x0 = -2
sig1[m] = 1
sig2[m + x0] = 1

freq = ft.fftfreq(n, centered=centered)
ftsig1 = ft.fft(sig1, centered=centered)
ftsig2 = ft.fft(sig2, centered=centered)
ftsig2_calc = ftsig1 * np.exp(-2 * np.pi * 1j * freq * x0)

image_product = ftsig1 * ftsig2.conj()
image_product_calc = np.exp(2 * np.pi * 1j * freq * x0) * np.abs(ftsig1) ** 2

image_product /= np.abs(image_product)
image_product_calc = np.exp(2 * np.pi * 1j * freq * x0)

cross_correlation = ft.ifft(image_product, centered=centered).real
cross_correlation_calc = np.zeros(n)
cross_correlation_calc[-x0] = 1

print(cross_correlation)
print(cross_correlation_calc)
print(-ft.fftfreq(n) * n)
