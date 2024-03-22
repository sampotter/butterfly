import colorcet as cc
import matplotlib.pyplot as plt; plt.ion()
import numpy as np
import pyvista as pv
import pyvistaqt as pvqt
import scipy.sparse

from matplotlib.colors import LogNorm

A = scipy.sparse.csr_matrix(
    (np.fromfile('A_data.bin', np.double),
     np.fromfile('A_indices.bin', np.intc),
     np.fromfile('A_indptr.bin', np.intc)))

M = scipy.sparse.csr_matrix(
    (np.fromfile('M_data.bin', np.double),
     np.fromfile('M_indices.bin', np.intc),
     np.fromfile('M_indptr.bin', np.intc)))

nodes = np.fromfile('nodes.bin', np.double).reshape(3, -1).T

mu = -0.001
U = scipy.sparse.linalg.eigsh(A, 10, M, mu, which='LM')[1]

poly_data = pv.PolyData(nodes)
poly_data['u'] = U[:, 9]
plotter = pvqt.BackgroundPlotter()
plotter.add_mesh(poly_data, cmap=cc.cm.gouldian, point_size=10)

x, y, z = nodes.T
phi = np.arccos(z)
theta = np.mod(np.arctan2(y, x), 2*np.pi)

plt.figure(figsize=(10, 6))
plt.scatter(theta, phi, 10, u, cmap=cc.cm.gouldian)
plt.xlim(0, 2*np.pi)
plt.ylim(np.pi, 0)
plt.gca().set_aspect('equal')
plt.tight_layout()
plt.xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi])
plt.yticks([np.pi, np.pi/2, 0])
plt.gca().set_xticklabels(['0', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'])
plt.gca().set_yticklabels(['0', r'$\pi/2$', r'$\pi$'][::-1])
plt.xlabel(r'$\theta$')
plt.ylabel(r'$\phi$')
plt.show()

# plt.figure(figsize=(11, 10))
# plt.title('A')
# plt.imshow(abs(A.A), cmap=cc.cm.gouldian, norm=LogNorm(clip=True))
# plt.colorbar()
# plt.gca().set_aspect('equal')
# plt.tight_layout()
# plt.show()

# plt.figure(figsize=(11, 10))
# plt.title('M')
# plt.imshow(abs(M.A), cmap=cc.cm.gouldian, norm=LogNorm(clip=True))
# plt.colorbar()
# plt.gca().set_aspect('equal')
# plt.tight_layout()
# plt.show()
