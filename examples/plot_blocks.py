#!/usr/bin/env python

import colorcet as cc
import matplotlib.pyplot as plt
import numpy as np
import sys

from matplotlib.patches import Rectangle

# plt.ion()

if len(sys.argv) != 2:
    print(f'usage: {sys.argv[0]} <blocks.txt>')
    exit()

f = open(sys.argv[1])
lines = f.readlines()
f.close()

Level, I0, I1, J0, J1, Type = np.array(
    [tuple(int(_) for _ in line.split()) for line in lines]
).T

assert(min(I0) == 0)
assert(min(J0) == 0)

max_level = Level.max()
default_bitmap_value = max_level + 1
num_rows = max(I1)
num_cols = max(J1)

# first make a bitmap---this is a useful quick plot and lets us make
# sure that we've hit every component of the matrix exactly once

bitmap = np.zeros((num_rows, num_cols), dtype=np.uint8)
bitmap[...] = default_bitmap_value
for i0, i1, j0, j1, level in zip(I0, I1, J0, J1, Level):
    assert (bitmap[i0:i1, j0:j1] == default_bitmap_value).all()
    bitmap[i0:i1, j0:j1] = level
assert (bitmap != default_bitmap_value).all()

plt.figure()
plt.imshow(bitmap, interpolation='none', cmap=cc.cm.rainbow)
plt.colorbar()
plt.show()

# now make a nicer block matrix plot

type2hatch = {
    7: '//', # BF_TYPE_MAT_DENSE_COMPLEX
    11: 'o',  # BF_MAT_PRODUCT (i.e.: butterfly factorization)
}

cmap = cc.cm.gouldian

fig, ax = plt.subplots()
for level, i0, i1, j0, j1, type_ in zip(Level, I0, I1, J0, J1, Type):
    di, dj = i1 - i0, j1 - j0
    fc = cmap(level/(max_level + 1))
    ec = cmap((level - 0.5)/(max_level + 1))
    if type_ not in type2hatch:
        print(type_)
    h = type2hatch[type_]
    rect = Rectangle((j0, i0), dj, di, linewidth=1, facecolor=fc,
                     edgecolor=ec, hatch=h)
    ax.add_patch(rect)
ax.set_aspect('equal')
ax.set_xlim(0, num_cols - 1)
ax.set_ylim(num_rows - 1, 0)
plt.show()
