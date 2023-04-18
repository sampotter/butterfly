import colorcet as cc
import matplotlib.pyplot as plt; plt.ion()
import numpy as np

from matplotlib.collections import EllipseCollection
from matplotlib.colors import LogNorm

xmin = -2
xmax = 2
ymin = -2
ymax = 2

nxEval = 401
nyEval = 401

ellipseData = np.fromfile('ellipseData.bin').reshape(-1, 5)

ellipseCenter = ellipseData[:, :2]
ellipseWidth = 2*ellipseData[:, 2]
ellipseHeight = 2*ellipseData[:, 3]
ellipseAngle = np.rad2deg(ellipseData[:, 4])

numEllipses = ellipseData.shape[0]

uIn = np.fromfile('uIn.bin', dtype=np.complex128).reshape(nxEval, nyEval)
uScatLu = np.fromfile('uScatLu.bin', dtype=np.complex128).reshape(nxEval, nyEval)
uLu = uIn + uScatLu

vmax = abs(np.real(uLu)).max()
vmin = -vmax

def myImshow(A, **kwargs):
    plt.imshow(np.rot90(A), extent=[xmin, xmax, ymin, ymax], **kwargs)

def plotEllipses():
    ellipseCollection = EllipseCollection(
        ellipseWidth, ellipseHeight,
        ellipseAngle,
        units='x',
        offsets=ellipseCenter,
        transOffset=plt.gca().transData,
        edgecolors='black',
        facecolors='white',
        linewidth=1,
    )
    plt.gca().add_collection(ellipseCollection)

cmap = cc.cm.gouldian

plt.figure(figsize=(13, 4))
plt.subplot(1, 3, 1)
myImshow(np.real(uIn), cmap=cmap, vmax=vmax, vmin=vmin)
plt.colorbar()
plotEllipses()
plt.subplot(1, 3, 2)
myImshow(np.real(uScatLu), cmap=cmap, vmax=vmax, vmin=vmin)
plt.colorbar()
plotEllipses()
plt.subplot(1, 3, 3)
myImshow(np.real(uLu), cmap=cmap, vmax=vmax, vmin=vmin)
plt.colorbar()
plotEllipses()
plt.tight_layout()
plt.show()

plt.figure(figsize=(12, 10))
myImshow(20*np.log10(abs(uLu)), cmap=cmap, vmax=0, vmin=-36)
plt.colorbar()
plotEllipses()
plt.tight_layout()
plt.show()
