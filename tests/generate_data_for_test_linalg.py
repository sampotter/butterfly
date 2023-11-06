#!/usr/bin/env python

'''Running this script the files required by test_linalg.c:

    sphere.obj
    sphere_Lam.txt
    sphere_Phi.txt

It doesn't need to be re-run each time the tests are compiled. The
purpose of including this file is to have it easily accessible in case
bugs are discovered and the test data needs to be tweaked.

'''

import numpy as np
import openmesh
import scipy.sparse
import trimesh

from array import array

SHOULD_PLOT = False
if SHOULD_PLOT:
    import colorcet as cc
    import pyvista as pv
    import pyvistaqt as pvqt
    import vtk

def get_fibonacci_sphere_points(n):
    x = np.linspace(1, -1, n)
    r, theta = np.sqrt(1 - x**2), np.pi*(np.sqrt(5) - 1)*np.arange(n)
    return np.array([x, r*np.cos(theta), r*np.sin(theta)]).T

def get_lbo_fem_discretization(mesh):
    L_data, M_data = array('d'), array('d')
    rows, cols = array('Q'), array('Q')

    for i, x in enumerate(mesh.vertices):
        # iterate over adjacent faces
        for j in (_ for _ in mesh.vertex_faces[i] if _ >= 0):
            # get the two face vertices opposite the current vertex
            i0, i1 = (_ for _ in mesh.faces[j] if _ != i)
            x0, x1 = mesh.vertices[i0], mesh.vertices[i1]

            # orthogonal projection of each vertex onto the opposite face edge
            d, d0, d1 = x1 - x0, x - x1, x0 - x
            y = x0 + d*np.dot(d, x - x0)/np.dot(d, d)
            y0 = x1 + d0*np.dot(d0, x0 - x1)/np.dot(d0, d0)
            y1 = x + d1*np.dot(d1, x1 - x)/np.dot(d1, d1)

            # compute gradients for hat functions centered at each vertex
            g = (x - y)/np.dot(x - y, x - y)
            g0 = (x0 - y0)/np.dot(x0 - y0, x0 - y0)
            g1 = (x1 - y1)/np.dot(x1 - y1, x1 - y1)

            # get triangle area
            area = np.linalg.norm(np.cross(x1 - x, x0 - x))/2

            # update L data
            L_data.extend([g.dot(g)*area, g.dot(g0)*area, g.dot(g1)*area])

            # update M data
            M_data.extend([area/6, area/12, area/12]) # TODO: double check this...

            # update COO index arrays
            rows.extend([i, i, i])
            cols.extend([i, i0, i1])

    L = scipy.sparse.coo_matrix((L_data, (rows, cols))).tocsr()
    M = scipy.sparse.coo_matrix((M_data, (rows, cols))).tocsr()

    # Both L and M should be numerically symmetric at this point. Go
    # and symmetrize them now...
    L = (L + L.T)/2
    M = (M + M.T)/2

    return L, M

if __name__ == '__main__':
    num_verts = 500
    lam_min = 50
    lam_max = 100

    verts = get_fibonacci_sphere_points(num_verts)
    faces = scipy.spatial.ConvexHull(verts).simplices
    mesh = trimesh.Trimesh(verts, faces)
    L, M = get_lbo_fem_discretization(mesh)

    Lam, Phi = scipy.sparse.linalg.eigsh(L, k=100, M=M, sigma=(lam_min + lam_max)/2)
    mask = np.logical_and(lam_min <= Lam, Lam <= lam_max)
    Lam = Lam[mask]
    Phi = Phi[:, mask]

    if SHOULD_PLOT:
        grid = pv.UnstructuredGrid({vtk.VTK_TRIANGLE: faces}, verts)
        grid['phi'] = Phi[:, 0]

        plotter = pvqt.BackgroundPlotter()
        plotter.add_mesh(grid, scalars='phi', cmap=cc.cm.coolwarm)

    with open('sphere.obj', 'w') as f:
        mesh.export(f, file_type='obj', digits=17, include_normals=True)

    np.savetxt('sphere_Lam.txt', Lam)
    np.savetxt('sphere_Phi.txt', Phi)
