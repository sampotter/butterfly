import colorcet as cc
import matplotlib.pyplot as plt
import numpy as np

f = open('blocks.txt')
lines = f.readlines()
f.close()

I0, I1, J0, J1, Level = np.array(
    [tuple(int(_) for _ in line.split()) for line in lines]
).T

max_level = Level.max()
default_bitmap_value = max_level + 1

bitmap = np.zeros((16384, 16384), dtype=np.uint8)
bitmap[...] = default_bitmap_value
for i0, i1, j0, j1, level in zip(I0, I1, J0, J1, Level):
    assert (bitmap[i0:i1, j0:j1] == default_bitmap_value).all()
    bitmap[i0:i1, j0:j1] = level

plt.figure()
plt.imshow(bitmap, interpolation='none', cmap=cc.cm.rainbow)
plt.colorbar()
plt.show()
