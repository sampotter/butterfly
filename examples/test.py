import colorcet as cc
import matplotlib.pyplot as plt
import numpy as np
import scipy.spatial
import scipy.special

from matplotlib.colors import hsv_to_rgb

np.random.seed(2345)

def Complex2HSV(z, rmin=None, rmax=None, hue_start=90):
    amp = np.abs(z)
    if rmin is None:
        rmin = amp.min()
    if rmax is None:
        rmax = amp.max()
    # get amplidude of z and limit to [rmin, rmax]
    amp = np.where(amp < rmin, rmin, amp)
    amp = np.where(amp > rmax, rmax, amp)
    ph = np.angle(z, deg=1) + hue_start
    # HSV are values in range [0,1]
    h = (ph % 360) / 360
    s = 0.85 * np.ones_like(h)
    v = (amp -rmin) / (rmax - rmin)
    return hsv_to_rgb(np.dstack((h,s,v)))

k = 3000

def evaluateKernel(arg):
    return (1j*scipy.special.j0(arg) - scipy.special.y0(arg))/4

def getKernelMatrix(X, Y):
    return evaluateKernel(k*scipy.spatial.distance_matrix(Y, X))

plt.ion()

first_srcCircPts = np.fromfile('first_srcCircPts.bin').reshape(-1, 2)
first_srcPts = np.fromfile('first_srcPts.bin').reshape(-1, 2)
first_tgtCircPts = np.fromfile('first_tgtCircPts.bin').reshape(-1, 2)

second_srcChildCircPts = np.fromfile('second_srcChildCircPts.bin').reshape(-1, 2)
second_tgtCircPts = np.fromfile('second_tgtCircPts.bin').reshape(-1, 2)
second_srcCircPts = np.fromfile('second_srcCircPts.bin').reshape(-1, 2)
second_tgtChildCircPts = np.fromfile('second_tgtChildCircPts.bin').reshape(-1, 2)

nx, ny = 401, 201
X, Y = np.meshgrid(np.linspace(-1, 1, nx), np.linspace(0, 1, ny), indexing='xy')
gridPts = np.array([X.ravel(), Y.ravel()]).T

########################################################################
# coefficients for charge distribution

V = np.random.randn(2) + 1j*np.random.randn(2)
v0 = np.random.randn() + 1j*np.random.randn()

########################################################################
# set up first problem

first_Q_or = v0 + first_srcPts@V

first_Q_eq = np.linalg.lstsq(
    getKernelMatrix(first_srcCircPts, first_tgtCircPts),
    getKernelMatrix(first_srcPts, first_tgtCircPts),
    rcond=None
)[0]@first_Q_or

first_Z_or = (getKernelMatrix(first_srcPts, gridPts)@first_Q_or).reshape(X.shape)
first_Z_eq = (getKernelMatrix(first_srcCircPts, gridPts)@first_Q_eq).reshape(X.shape)

first_errors = -np.log10(np.maximum(1e-16, abs(first_Z_or - first_Z_eq)))

########################################################################
# set up second problem

second_Q_or = v0 + second_srcChildCircPts@V

second_Q_eq = np.linalg.lstsq(
    getKernelMatrix(second_srcCircPts, second_tgtChildCircPts),
    getKernelMatrix(second_srcChildCircPts, second_tgtChildCircPts),
    rcond=None
)[0]@second_Q_or

second_Z_or = (getKernelMatrix(second_srcChildCircPts, gridPts)@second_Q_or).reshape(X.shape)
second_Z_eq = (getKernelMatrix(second_srcCircPts, gridPts)@second_Q_eq).reshape(X.shape)

second_errors = -np.log10(np.maximum(1e-16, abs(second_Z_or - second_Z_eq)))

########################################################################
# make plots

extent = [-1, 1, 1, 0]

plt.figure(figsize=(16, 9))

# first problem (original)

plt.subplot(2, 3, 1)
plt.imshow(Complex2HSV(first_Z_or, rmin=0), extent=extent, zorder=1)
plt.scatter(*first_srcPts.T, s=0.5, c='r', zorder=2)
plt.scatter(*first_srcCircPts.T, s=0.5, c='pink', zorder=3)
plt.scatter(*first_tgtCircPts.T, s=0.5, c='b', zorder=4)
plt.xlim(-1, 1)
plt.ylim(0, 1)
plt.gca().set_aspect('equal')
plt.xlabel(r'$x$')
plt.xlabel(r'$y$')
plt.title('first problem (original field)')

# first problem (equivalent---after shift)

plt.subplot(2, 3, 2)
plt.imshow(Complex2HSV(first_Z_eq, rmin=0), extent=extent, zorder=1)
plt.scatter(*first_srcPts.T, s=0.5, c='r', zorder=2)
plt.scatter(*first_srcCircPts.T, s=0.5, c='pink', zorder=3)
plt.scatter(*first_tgtCircPts.T, s=0.5, c='b', zorder=4)
plt.xlim(-1, 1)
plt.ylim(0, 1)
plt.gca().set_aspect('equal')
plt.xlabel(r'$x$')
plt.xlabel(r'$y$')
plt.title('first problem (shifted field)')

# first problem (errors)

plt.subplot(2, 3, 3)
plt.imshow(first_errors, cmap=cc.cm.rainbow, extent=extent, zorder=1)
plt.colorbar()
plt.scatter(*first_srcPts.T, s=0.5, c='r', zorder=2)
plt.scatter(*first_srcCircPts.T, s=0.5, c='pink', zorder=3)
plt.scatter(*first_tgtCircPts.T, s=0.5, c='b', zorder=4)
plt.xlim(-1, 1)
plt.ylim(0, 1)
plt.gca().set_aspect('equal')
plt.xlabel(r'$x$')
plt.xlabel(r'$y$')
plt.title('first problem (shifted field)')

# second problem (original)

plt.subplot(2, 3, 4)
plt.imshow(Complex2HSV(second_Z_or, rmin=0), extent=extent, zorder=1)
plt.scatter(*second_srcChildCircPts.T, s=0.5, c='r', zorder=2)
plt.scatter(*second_srcCircPts.T, s=0.5, c='pink', zorder=3)
plt.scatter(*second_tgtCircPts.T, s=0.5, c='b', zorder=4)
plt.scatter(*second_tgtChildCircPts.T, s=0.5, c='cyan', zorder=5)
plt.xlim(-1, 1)
plt.ylim(0, 1)
plt.gca().set_aspect('equal')
plt.xlabel(r'$x$')
plt.xlabel(r'$y$')
plt.title('second problem (shifted field)')

# second problem (equivalent---after shift)

plt.subplot(2, 3, 5)
plt.imshow(Complex2HSV(second_Z_eq, rmin=0), extent=extent, zorder=1)
plt.scatter(*second_srcChildCircPts.T, s=0.5, c='r', zorder=2)
plt.scatter(*second_srcCircPts.T, s=0.5, c='pink', zorder=3)
plt.scatter(*second_tgtCircPts.T, s=0.5, c='b', zorder=4)
plt.scatter(*second_tgtChildCircPts.T, s=0.5, c='cyan', zorder=5)
plt.xlim(-1, 1)
plt.ylim(0, 1)
plt.gca().set_aspect('equal')
plt.xlabel(r'$x$')
plt.xlabel(r'$y$')
plt.title('second problem (shifted field)')

# second problem (errors)

plt.subplot(2, 3, 6)
plt.imshow(second_errors, cmap=cc.cm.rainbow, extent=extent, zorder=1)
plt.colorbar()
plt.scatter(*second_srcChildCircPts.T, s=0.5, c='r', zorder=2)
plt.scatter(*second_srcCircPts.T, s=0.5, c='pink', zorder=3)
plt.scatter(*second_tgtCircPts.T, s=0.5, c='b', zorder=4)
plt.scatter(*second_tgtChildCircPts.T, s=0.5, c='cyan', zorder=5)
plt.xlim(-1, 1)
plt.ylim(0, 1)
plt.gca().set_aspect('equal')
plt.xlabel(r'$x$')
plt.xlabel(r'$y$')
plt.title('first problem (shifted field)')

plt.tight_layout()
