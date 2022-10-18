import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse
import scipy.sparse.linalg

plt.ion()

from PIL import Image
from util import complex_to_hsv

n = 2048

rowInd = np.fromfile('rowInd.bin', dtype=np.uintp)
colInd = np.fromfile('colInd.bin', dtype=np.uintp)
value = np.fromfile('value.bin', dtype=np.complex128)

A_lift = scipy.sparse.coo_matrix((value, (rowInd, colInd)))

# vmax = abs(A_lift).max()
# vmin = -vmax

# A_lift_hsv = complex_to_hsv(A_lift.toarray())

# A_lift_hsv_img = Image.fromarray(A_lift_hsv, 'RGB')

# plt.figure()
# plt.imshow()
# plt.show()

plt.figure()
plt.spy(A_lift, marker='.', markersize=0.1, c='k')
plt.show()

lu = scipy.sparse.linalg.splu(A_lift.tocsc())

print(f'n = {n}')
print(f'A.size = {n**2}')
print(f'A_BF.size = {value.size}')
print(f'nnz(L) = {lu.L.nnz}')
print(f'nnz(U) = {lu.U.nnz}')

plt.figure()
plt.scatter(np.arange(lu.perm_r.size), lu.perm_r, marker='.', s=1, c='k')
plt.gca().set_aspect('equal')
plt.show()
