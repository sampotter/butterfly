#!/usr/bin/env python

import colorcet as cc
import matplotlib.pyplot as plt
import numpy as np

from util import complex_to_hsv

plt.ion()

A_dense = np.fromfile('A_dense.bin', dtype=np.complex128)

num_points = np.sqrt(A_dense.size)
assert np.modf(np.sqrt(A_dense.size))[0] == 0
num_points = int(num_points)

A_dense = A_dense.reshape(num_points, num_points)

A_BF_dense = np.fromfile('A_BF_dense.bin', dtype=np.complex128)
A_BF_dense = A_BF_dense.reshape(num_points, num_points)

error = A_dense - A_BF_dense

log_error = -np.log10(np.maximum(1e-16, abs(error)))

print(f'minimum log10 error: {log_error.min()}')

plt.figure()
plt.imshow(complex_to_hsv(A_BF_dense), interpolation='none')
plt.title('BF')
plt.show()

plt.figure()
plt.imshow(complex_to_hsv(A_dense), interpolation='none')
plt.title('Dense')
plt.show()

plt.figure()
plt.imshow(log_error, cmap=cc.cm.rainbow, interpolation='none')
plt.colorbar()
plt.show(block=True)
