import colorcet as cc
import matplotlib.pyplot as plt
import numpy as np

from matplotlib.colors import LogNorm

plt.figure()
plt.plot(abs(np.fromfile('uDenseCheck.bin', dtype=np.complex128)))
plt.show()

xmin = -1
xmax = 1
ymin = -1
ymax = 1

nxEval = 401
nyEval = 401

uIn = np.fromfile('uIn.bin', dtype=np.complex128).reshape(nxEval, nyEval)
uScatDense = np.fromfile('uScatDense.bin', dtype=np.complex128).reshape(nxEval, nyEval)
uDense = uIn + uScatDense

vmax = abs(np.real(uDense)).max()
vmin = -vmax

def myImshow(A, **kwargs):
    plt.imshow(np.rot90(A), extent=[xmin, xmax, ymin, ymax], **kwargs)

plt.figure(figsize=(13, 4))
plt.subplot(1, 3, 1)
myImshow(np.real(uIn), cmap=cc.cm.coolwarm, vmax=vmax, vmin=vmin)
plt.colorbar()
plt.subplot(1, 3, 2)
myImshow(np.real(uScatDense), cmap=cc.cm.coolwarm, vmax=vmax, vmin=vmin)
plt.colorbar()
plt.subplot(1, 3, 3)
myImshow(np.real(uDense), cmap=cc.cm.coolwarm, vmax=vmax, vmin=vmin)
plt.colorbar()
plt.tight_layout()
plt.show()

kwargs = {
    'cmap': cc.cm.coolwarm,
    'norm': LogNorm(vmin=1e-5, vmax=1)
}
plt.figure(figsize=(9, 8))
myImshow(abs(uDense), **kwargs)
plt.colorbar()
plt.tight_layout()
plt.show()
