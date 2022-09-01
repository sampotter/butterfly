import colorcet as cc
import itertools as it
import matplotlib.pyplot as plt
import numpy as np
import scipy.special

from util import complex_to_hsv as Complex2HSV

tgt_pts = np.fromfile('tgt_pts.bin', dtype=np.float64).reshape(-1, 2)
src_pts = np.fromfile('src_pts.bin', dtype=np.float64).reshape(-1, 2)

tgt_circ_pts = np.fromfile('tgt_circ_pts.bin', dtype=np.float64).reshape(-1, 2)
src_circ_pts = np.fromfile('src_circ_pts.bin', dtype=np.float64).reshape(-1, 2)

m, n, p = tgt_pts.shape[0], src_pts.shape[0], src_circ_pts.shape[0]
assert p == tgt_circ_pts.shape[0]

plt.figure()
plt.scatter(*tgt_pts.T, s=3, c=np.linspace(0, 0.5, m), marker='.', cmap=cc.cm.bmw)
plt.scatter(*src_pts.T, s=3, c=np.linspace(0.5, 1, m), marker='.', cmap=cc.cm.bmw)
plt.colorbar()
plt.scatter(*tgt_circ_pts.T, s=2, c='k', marker='.')
plt.scatter(*src_circ_pts.T, s=2, c='k', marker='.')
plt.gca().set_aspect('equal')
plt.show()

Z_gt = np.fromfile('Z_gt.bin', dtype=np.complex128).reshape(m, n)

U_gt = np.fromfile('U_gt.bin', dtype=np.complex128).reshape(m, m)
S_gt = np.fromfile('S_gt.bin', dtype=np.float64).reshape(m)
Vt_gt = np.fromfile('Vt_gt.bin', dtype=np.complex128).reshape(n, n)

Z1 = np.fromfile('Z1.bin', dtype=np.complex128).reshape(p, n)
Z2 = np.fromfile('Z2.bin', dtype=np.complex128).reshape(p, p)
Z3 = np.fromfile('Z3.bin', dtype=np.complex128).reshape(m, p)

U2 = np.fromfile('U2.bin', dtype=np.complex128).reshape(p, p)
S2 = np.fromfile('S2.bin', dtype=np.float64).reshape(p)
Vt2 = np.fromfile('Vt2.bin', dtype=np.complex128).reshape(p, p)

x = np.fromfile('x.bin', dtype=np.complex128).reshape(n, -1)
b = np.fromfile('b.bin', dtype=np.complex128).reshape(m, -1)
b_gt = np.fromfile('b_gt.bin', dtype=np.complex128).reshape(m, -1)

num_trials = x.shape[1]
assert b.shape[1] == num_trials and b_gt.shape[1] == num_trials

############################################################################
# quick test for computing shift matrices

Z_or = np.fromfile('Z_or.bin', dtype=np.complex128).reshape(17, 17)
Z_eq = np.fromfile('Z_eq.bin', dtype=np.complex128).reshape(17, 17)
# Z_eq_pinv = np.fromfile('Z_eq_pinv.bin', dtype=np.complex128).reshape(19, 19)
# Z_eq_pinv_gt = np.linalg.pinv(Z_eq)
Z_shift = np.fromfile('Z_shift.bin', dtype=np.complex128).reshape(17, 17)
Z_shift_gt = np.linalg.lstsq(Z_eq, Z_or, rcond=None)[0]

plt.figure(figsize=(6, 8))
plt.subplot(2, 2, 1)
plt.imshow(Complex2HSV(Z_or, rmin=0))
plt.subplot(2, 2, 2)
plt.imshow(Complex2HSV(Z_eq, rmin=0))
plt.subplot(2, 2, 3)
plt.imshow(Complex2HSV(Z_shift, rmin=0))
plt.subplot(2, 2, 4)
plt.imshow(Complex2HSV(Z_shift_gt, rmin=0))
plt.tight_layout()
plt.show()

############################################################################
# check sparsity pattern of butterfly factors

import glob

paths = glob.glob("ij*.txt")

IJs = dict()

for path in paths:
    level = int(path.split('.')[0][2:])
    IJs[level] = np.loadtxt(path).astype(int)

mats = dict()
for level, IJ in IJs.items():
    m, n = IJ.max(0) + 1
    mat = np.empty((m, n), dtype=np.float64)
    mat[...] = np.nan
    for c, (i, j) in enumerate(IJ):
        mat[i, j] = c
    mats[level] = mat

plt.figure()
for level, mat in mats.items():
    plt.subplot(1, len(mats), len(mats) - level + 1)
    plt.imshow(mat, cmap=cc.cm.colorwheel)
    plt.gca().set_aspect('equal')
    plt.title('level = %d' % (level,))
plt.show()
