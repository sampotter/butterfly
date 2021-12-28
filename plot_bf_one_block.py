#!/usr/bin/env python

import colorcet as cc
import matplotlib.pyplot as plt
import numpy as np
import scipy.spatial
import scipy.special

from glob import glob
from util import complex_to_hsv
from pathlib import Path

SAVE_PLOTS = True

if SAVE_PLOTS:
    plots_path = Path('plots')
    plots_path.mkdir()

f = open('info.txt', 'r')
info = dict()
for line in f.readlines():
    key, value = line.split()
    info[key] = value
f.close()

info['numSrcPts'] = int(info['numSrcPts'])
info['numTgtPts'] = int(info['numTgtPts'])
info['K'] = float(info['K'])

for key, value in info.items():
    print(f'{key}: {value}')

K = info['K'] # wavenumber

M, N = info['numTgtPts'], info['numSrcPts']

# load groundtruth kernel matrix
Z_gt = np.fromfile('Z_gt.bin', dtype=np.complex128).reshape(M, N)

# load source and target points
srcPts = np.fromfile('srcPts.bin', dtype=np.float64).reshape(-1, 2)
tgtPts = np.fromfile('tgtPts.bin', dtype=np.float64).reshape(-1, 2)

sortedSrcPts = srcPts[np.argsort(srcPts[:, 1])]
sortedTgtPts = tgtPts[np.argsort(tgtPts[:, 1])]

########################################################################
# set up test problem

np.random.seed(1234)

# q = np.random.randn(3) + 1j*np.random.randn(3)
# Q = q[0] + srcPts@q[1:]

Q = np.random.randn(M) + 1j*np.random.randn(M)

def getZ(srcPts, tgtPts):
    arg = K*scipy.spatial.distance_matrix(tgtPts, srcPts)
    return (1j*scipy.special.j0(arg) - scipy.special.y0(arg))/4

def getshift(srcPtsOrig, srcPtsEquiv, tgtPts):
    Z_or = getZ(srcPtsOrig, tgtPts)
    Z_eq = getZ(srcPtsEquiv, tgtPts)
    return np.linalg.lstsq(Z_eq, Z_or, rcond=None)[0]

def getphi(srcPts, tgtPts, Q):
    arg = K*scipy.spatial.distance_matrix(tgtPts, srcPts)
    try:
        phi = (1j*scipy.special.j0(arg) - scipy.special.y0(arg))@Q/4
    except:
        import ipdb; ipdb.set_trace()
    return phi

########################################################################
# function for plotting the difference between original and
# re-expanded field

def plotdiff(srcPtsOrig, srcPtsEquiv, tgtPts, qOrig, block, h=0.0025, bbox=None):
    if bbox is None:
        xmin, ymin = np.array([
            srcPtsOrig.min(0), srcPtsEquiv.min(0), tgtPts.min(0)
        ]).min(0)

        xmax, ymax = np.array([
            srcPtsOrig.max(0), srcPtsEquiv.max(0), tgtPts.max(0)
        ]).max(0)

        xmin -= 0.1
        xmax += 0.1
        ymin -= 0.1
        ymax += 0.1
    else:
        xmin, xmax, ymin, ymax = bbox

    dx, dy = xmax - xmin, ymax - ymin
    extent = [xmin, xmax, ymax, ymin]
    nx, ny = int(np.ceil(dx/h)), int(np.ceil(dy/h))
    X, Y = np.meshgrid(
        np.linspace(xmin, xmax, nx),
        np.linspace(ymin, ymax, ny),
        indexing='xy'
    )
    srcPtsEval = np.array([X.ravel(), Y.ravel()]).T

    PhiOrig = getphi(srcPtsOrig, srcPtsEval, qOrig).reshape(X.shape)

    qEquiv = block@qOrig
    PhiEquiv = getphi(srcPtsEquiv, srcPtsEval, qEquiv).reshape(X.shape)

    rel_diff = (PhiEquiv - PhiOrig)/max(abs(PhiEquiv).max(), abs(PhiOrig).max())
    digits = -np.log10(np.maximum(1e-16, abs(rel_diff)))

    fig = plt.figure(figsize=(10, 8))

    ax = fig.add_subplot(3, 1, 1)
    ax.imshow(complex_to_hsv(PhiOrig), extent=extent)
    ax.scatter(*srcPtsOrig.T, s=1, c='cyan', zorder=1)
    ax.scatter(*srcPtsEquiv.T, s=1, c='orange', zorder=2)
    ax.scatter(*tgtPts.T, s=1, c='pink', zorder=3)
    ax.plot(*sortedSrcPts.T, linewidth=1, linestyle='--', c='white', zorder=4)
    ax.plot(*sortedTgtPts.T, linewidth=1, linestyle='--', c='white', zorder=4)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect('equal')

    ax = fig.add_subplot(3, 1, 2)
    ax.imshow(complex_to_hsv(PhiEquiv), extent=extent)
    ax.scatter(*srcPtsOrig.T, s=1, c='cyan', zorder=1)
    ax.scatter(*srcPtsEquiv.T, s=1, c='orange', zorder=2)
    ax.scatter(*tgtPts.T, s=1, c='pink', zorder=3)
    ax.plot(*sortedSrcPts.T, linewidth=1, linestyle='--', c='white', zorder=4)
    ax.plot(*sortedTgtPts.T, linewidth=1, linestyle='--', c='white', zorder=4)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect('equal')

    ax = fig.add_subplot(3, 1, 3)
    im = ax.imshow(digits, extent=extent, cmap=cc.cm.rainbow, vmin=0,
                   vmax=16, zorder=0)
    fig.colorbar(im, ax=ax)
    ax.scatter(*srcPtsOrig.T, s=1, c='cyan', zorder=1)
    ax.scatter(*srcPtsEquiv.T, s=1, c='orange', zorder=2)
    ax.scatter(*tgtPts.T, s=1, c='pink', zorder=3)
    ax.plot(*sortedSrcPts.T, linewidth=1, linestyle='--', c='white', zorder=4)
    ax.plot(*sortedTgtPts.T, linewidth=1, linestyle='--', c='white', zorder=4)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect('equal')

    fig.tight_layout()

    return fig

########################################################################
# load factors

xmin, ymin = np.minimum(srcPts.min(0), tgtPts.min(0)) - 0.25
xmax, ymax = np.maximum(srcPts.max(0), tgtPts.max(0)) + 0.25
bbox = (xmin, xmax, ymin, ymax)

nnz = 0 # count number of nonzero entries in the factor blocks

Qs = [Q.copy()]

factors = []

frame = 0

paths = glob('factor*')
for factorNum, path in enumerate(paths):
    print(f'loading {path}')

    Qprev = Qs[-1]

    f = open(f'{path}/info.txt')
    factor_info = dict()
    for line in f.readlines():
        key, value = line.split()
        factor_info[key] = value
    f.close()

    for key in ['numBlockRows', 'numBlockCols', 'numBlocks']:
        factor_info[key] = int(factor_info[key])

    rowInd = np.fromfile(f'{path}/rowInd.bin', dtype=np.uintp)
    colInd = np.fromfile(f'{path}/colInd.bin', dtype=np.uintp)

    rowOffset = np.fromfile(f'{path}/rowOffset.bin', dtype=np.uintp)
    colOffset = np.fromfile(f'{path}/colOffset.bin', dtype=np.uintp)

    factor = np.zeros((rowOffset[-1], colOffset[-1]), dtype=np.complex128)

    for k, (i, j) in enumerate(zip(rowInd, colInd)):
        print(f'* block #{k}')

        i0, i1 = rowOffset[i], rowOffset[int(i + 1)]
        j0, j1 = colOffset[j], colOffset[int(j + 1)]

        assert (factor[i0:i1, j0:j1] == 0).all()

        block = np.fromfile(f'{path}/block{k}.bin', dtype=np.complex128)
        block = block.reshape(i1 - i0, j1 - j0)

        if path != paths[-1]:
            srcPtsOrig = np.fromfile(f'{path}/srcPtsOrig{k}.bin', dtype=np.float64)
            srcPtsOrig = srcPtsOrig.reshape(-1, 2)

            srcPtsEquiv = np.fromfile(f'{path}/srcPtsEquiv{k}.bin', dtype=np.float64)
            srcPtsEquiv = srcPtsEquiv.reshape(-1, 2)

            tgtPts = np.fromfile(f'{path}/tgtPts{k}.bin', dtype=np.float64)
            tgtPts = tgtPts.reshape(-1, 2)

            fig = plotdiff(srcPtsOrig, srcPtsEquiv, tgtPts, Qprev[j0:j1], block, bbox=bbox)
            if SAVE_PLOTS:
                img_path = plots_path/f'frame{frame}_factor{factorNum}_block{k}_i{i}_j{j}.png'
                fig.savefig(img_path)
            else:
                fig.show()

        nnz += block.size

        factor[i0:i1, j0:j1] = block

        frame += 1

    factors.append(factor)
    Qs.append(factor@Qprev)

Z_butterfly = factors[0]
for i in range(1, len(factors)):
    Z_butterfly = factors[i]@Z_butterfly

########################################################################
# stats

print(f'Z_gt elts: {Z_gt.size}')
print(f'factor elts: {nnz}')

########################################################################
# make plots

# plot groundtruth kernel matrix and product of butterfly factors

plt.figure(figsize=(11, 4))
plt.subplot(1, 3, 1)
plt.imshow(complex_to_hsv(Z_gt, rmin=0))
plt.gca().set_aspect('equal')
plt.title('Z (groundtruth)')
plt.subplot(1, 3, 2)
plt.imshow(complex_to_hsv(Z_butterfly, rmin=0))
plt.gca().set_aspect('equal')
plt.title('Z (butterfly)')
plt.tight_layout()
plt.subplot(1, 3, 3)
plt.imshow(-np.log10(np.maximum(np.finfo(np.float64).eps,
                                abs(Z_butterfly - Z_gt)/abs(Z_gt))),
           cmap=cc.cm.rainbow)
plt.colorbar()
plt.gca().set_aspect('equal')
if SAVE_PLOTS:
    plt.savefig(plots_path/'compare.png')
else:
    plt.show()

# plot blocks

for k, Z in enumerate(factors):
    plt.figure()
    plt.imshow(complex_to_hsv(Z, rmin=0))
    plt.gca().set_aspect('equal')
    plt.title(f'Factor #{k}')
    if SAVE_PLOTS:
        plt.savefig(plots_path/f'factor{k}.png')
    else:
        plt.show()
