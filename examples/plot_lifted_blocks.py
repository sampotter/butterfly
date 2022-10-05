#!/usr/bin/env python

import colorcet as cc
import matplotlib.pyplot as plt
import numpy as np
import sys

from matplotlib.patches import Rectangle

if len(sys.argv) != 2:
    print(f'usage: {sys.argv[0]} <blocks.txt>')
    sys.exit()

f = open(sys.argv[1])
lines = f.readlines()
f.close()

I0, J0, M, N, Type = np.array(
    [tuple(int(_) for _ in line.split()) for line in lines]
).T

I1, J1 = I0 + M, J0 + N
num_rows, num_cols = I1.max(), J1.max()

min_t = Type.min()
max_t = Type.max()

cmap = cc.cm.glasbey

MAT = 0
MAT_BLOCK = 1
MAT_BLOCK_COO = 2
MAT_BLOCK_DENSE = 3
MAT_BLOCK_DIAG = 4
MAT_COO_COMPLEX = 5
MAT_COO_REAL = 6
MAT_DENSE_COMPLEX = 7
MAT_DENSE_REAL = 8
MAT_DIAG_REAL = 9
MAT_GIVENS_COMPLEX = 10
MAT_PRODUCT = 11
MAT_SUM = 12
MAT_ZERO = 13

plt.figure()
rects = []
for i0, j0, m, n, t in zip(I0, J0, M, N, Type):
    if t == MAT_ZERO:
        continue
    fc = cmap((t - min_t)/(max_t - min_t))
    rect = Rectangle((j0, i0), n, m, linewidth=1, facecolor=fc, edgecolor='k')
    plt.gca().add_patch(rect)
plt.gca().set_aspect('equal')
plt.xlim(0, num_cols - 1)
plt.ylim(0, num_rows - 1)
plt.gca().invert_yaxis()
plt.show()
