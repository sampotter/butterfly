import colorcet as cc
import itertools as it
import matplotlib.pyplot as plt
import numpy as np
import scipy.special

from matplotlib.colors import hsv_to_rgb

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

def Complex2HSV(z, rmin, rmax, hue_start=90):
    # get amplidude of z and limit to [rmin, rmax]
    amp = np.abs(z)
    amp = np.where(amp < rmin, rmin, amp)
    amp = np.where(amp > rmax, rmax, amp)
    ph = np.angle(z, deg=1) + hue_start
    # HSV are values in range [0,1]
    h = (ph % 360) / 360
    s = 0.85 * np.ones_like(h)
    v = (amp -rmin) / (rmax - rmin)
    return hsv_to_rgb(np.dstack((h,s,v)))


Z_or = np.fromfile('Z_or.bin', dtype=np.complex128).reshape(19, 34)
Z_eq = np.fromfile('Z_eq.bin', dtype=np.complex128).reshape(19, 19)
Z_eq_pinv = np.fromfile('Z_eq_pinv.bin', dtype=np.complex128).reshape(19, 19)

mask = S > 1e-10
_ = (VH[mask].conj().T@(U[:, mask]/S[mask]).conj().T)@Z_or
plt.imshow(Complex2HSV(_, 0, abs(_).max()))
