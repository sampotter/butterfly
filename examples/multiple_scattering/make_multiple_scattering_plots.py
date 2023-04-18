#!/usr/bin/env python

import argparse
import colorcet as cc
import matplotlib.pyplot as plt; plt.ion()
import numpy as np

from matplotlib.collections import EllipseCollection
from matplotlib.colors import LogNorm

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('uIn', type=str)
    parser.add_argument('uScat', type=str)
    parser.add_argument('--xmin', type=float)
    parser.add_argument('--xmax', type=float)
    parser.add_argument('--ymin', type=float)
    parser.add_argument('--ymax', type=float)
    parser.add_argument('--nx', type=int)
    parser.add_argument('--ny', type=int)
    args = parser.parse_args()

    ellipseData = np.fromfile('ellipseData.bin').reshape(-1, 5)

    ellipseCenter = ellipseData[:, :2]
    ellipseWidth = 2*ellipseData[:, 2]
    ellipseHeight = 2*ellipseData[:, 3]
    ellipseAngle = np.rad2deg(ellipseData[:, 4])

    numEllipses = ellipseData.shape[0]

    uIn = np.fromfile(args.uIn, dtype=np.complex128).reshape(args.nx, args.ny)
    uScat = np.fromfile(args.uScat, dtype=np.complex128).reshape(args.nx, args.ny)
    u = uIn + uScat

    vmax = abs(np.real(u)).max()
    vmin = -vmax

    def myImshow(A, **kwargs):
        plt.imshow(np.rot90(A), extent=[args.xmin, args.xmax, args.ymin, args.ymax], **kwargs)

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
    myImshow(np.real(uScat), cmap=cmap, vmax=vmax, vmin=vmin)
    plt.colorbar()
    plotEllipses()
    plt.subplot(1, 3, 3)
    myImshow(np.real(u), cmap=cmap, vmax=vmax, vmin=vmin)
    plt.colorbar()
    plotEllipses()
    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(12, 10))
    myImshow(20*np.log10(abs(u)), cmap=cmap, vmax=0, vmin=-36)
    plt.colorbar()
    plotEllipses()
    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(12, 10))
    myImshow(np.real(u), cmap=cmap, vmin=vmin, vmax=vmax)
    plt.colorbar()
    plotEllipses()
    plt.tight_layout()
    plt.show()
