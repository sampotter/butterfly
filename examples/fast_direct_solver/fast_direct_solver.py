#!/usr/bin/env python3

'''OVERVIEW
--------

This script contains a prototype of a hierarchical LU decomposition
based on the butterfly factorization.

What's novel about it? Let's see:

- ...

NOTES & QUESTIONS
-----------------

Should I be butterflying A21*inv(A11) and inv(A11)*A12 instead of
A21*inv(A11)*A12? -> YES. Can't get a good complexity otherwise. But a
problem with this is that it's not clear how to prove that these
matrices are butterfliable (although they appear to be?).

I'm not sure if the pseudoinverses of A21 and A12 are useful here. But
it might be possible to construct an "pseudoinverse butterfly" by
trying to run the physics backwards. Can we get the original field
back? This is tantamount to solving a ton of small inverse
problems... not sure.

For handy reference, here's a derivation of the inversion formula for
an invertible 2x2 block matrix:

Start with:

   A x = b

=> [ Id           0  ] [ A11 0                      ] [ Id inv(A11)*A12 ] [ x1 ] = [ y1 ]
   [ A21*inv(A11) Id ] [ 0   A22 - A21*inv(A11)*A12 ] [ 0  Id           ] [ x2 ] = [ y2 ]

=> [ x1 ] = [ Id -inv(A11)*A12 ] [ inv(A11) 0                           ] [  Id              ] [ y1 ]
   [ x2 ] = [     Id           ] [ 0        inv(A22 - A21*inv(A11)*A12) ] [ -A21*inv(A11) Id ] [ y2 ]

=> [ x1 ] = [ Id -inv(A11)*A12 ] [ inv(A11) 0                           ] [  y1                   ]
   [ x2 ] = [     Id           ] [ 0        inv(A22 - A21*inv(A11)*A12) ] [ -A21*inv(A11)*y1 + y2 ]

=> [ x1 ] = [ Id -inv(A11)*A12 ] [ inv(A11)*y1                                        ]
   [ x2 ] = [     Id           ] [ inv(A22 - A21*inv(A11)*A12)*(y2 - A21*inv(A11)*y1) ]

=> [ x1 ] = [ inv(A11)*y1 - inv(A11)*A12*inv(A22 - A21*inv(A11)*A12)*(y2 - A21*inv(A11)*y1) ]
   [ x2 ] = [ inv(A22 - A21*inv(A11)*A12)*(y2 - A21*inv(A11)*y1)                            ]

Should verify this using some random invertible matrix to make sure
it's bug-free...

'''

MAKE_PLOTS = False

import sys
sys.path.insert(-1, '..') # for util
sys.path.insert(-1, '../../wrappers/python') # for butterfly

import argparse
import numpy as np
import scipy.linalg
import scipy.sparse

from sklearn.utils.extmath import randomized_svd

import butterfly as bf

from util import complex_to_hsv

if MAKE_PLOTS:
    import matplotlib.pyplot as plt

    from matplotlib.collections import PatchCollection
    from matplotlib.patches import Rectangle

    def run_from_ipython():
        '''Returns True if this script is being run from IPython.'''
        try:
            __IPYTHON__
            return True
        except:
            return False

    # If we're running this script from IPython, turn matplotlib's
    # interactive mode on so plots don't block (if we're running the
    # script from the command line, we'd like the plots to block so we
    # have a chance to see them before they disappear...)
    if run_from_ipython(): plt.ion()

# simulation parameters
seed = 0
xmin, xmax = -1, 1
ymin, ymax = -1, 1
xsrc = bf.Points2.from_point((-0.5, 2))
k = 50
r = 0.2 # spacing between ellipses
axlim = (0.02, 0.1) # range of ellipse semi axes
h = 0.01 # point spacing
KR_order = 6 # Kapur-Rokhlin quadrature order

# seed bf's PRNG
bf.seed(seed)

# set up 2D Helmholtz kernel
helm = bf.Helm2(k, bf.LayerPot.Sp)

# use Poisson disk sampling to generate ellipse centers
bbox = bf.Bbox2(xmin + r, xmax - r, ymin + r, ymax - r)
centers = bf.Points2.sample_poisson_disk(bbox, r)
num_ellipses = len(centers)

# randomly sample ellipses
ellipses = []
for i in range(num_ellipses):
    a, b = bf.sample_uniform(*axlim, n=2)
    theta = bf.sample_uniform(0, 2*np.pi)
    ellipses.append(bf.Ellipse(max(a, b), min(a, b), centers[i], theta))
print(f'set up {len(ellipses)} random ellipses')

# discretize the ellipses to get the problem geometry
offsets = [0]
X, W, N = bf.Points2(), bf.RealArray(), bf.Vectors2()
for ellipse in ellipses:
    n = int(np.ceil(ellipse.perimeter/h))
    offsets.append(offsets[-1] + n)
    X_, _, N_, W_ = ellipse.sample_linspaced(n)
    X.extend(X_)
    N.extend(N_)
    W.extend(W_)
n = len(X)
print(f'set up discretization with {n} points')

# build quadtree on points w/ normals
quadtree = bf.Quadtree.from_points_and_normals(X, N)

print(np.array(W))

# get permuted copy of quadrature weights
rev_perm = quadtree.perm.get_reverse()
W_perm = W.copy()
print(np.array(W_perm))
W_perm.permute(rev_perm)

print(np.array(W_perm))

# compute the incident field
phi_in = helm.get_kernel_matrix(xsrc, Xtgt=X, Ntgt=N)
A = bf.FacHelm2.make_multilevel(helm, quadtree)

# set up multilevel butterfly factorization
helm.apply_block_KR_correction(A, offsets, KR_order, X, N, tree=quadtree)
A.scale_cols(W_perm)
A += bf.MatIdentity(n)/2

def get_block_inds_for_split(nodes):
    '''Find the dimension we should split along by fitting a tight
    bounding box to `nodes`' points and checking which dimension is
    longer. Then partition `nodes` into two groups depending on which
    half of the bounding box their centroids lie, giving two index
    vectors---`block_inds[0]` and `block_inds[1]` partitioning
    `range(len(nodes))`. The index vectors are returned so that
    `len(block_inds[0]) <= len(block_inds[1])`.

    '''
    # Find the "long" dimension:
    X = bf.Points2()
    for node in nodes:
        X.extend(node.get_points())
    long_dim = np.ptp(X, axis=0).argmax()

    # Get the points' centroid, which we'll use to partition the
    # nodes, comparing it with the quadtree centroids along the long
    # dimension.
    centroid = np.mean(X, axis=0)

    # Find the index vectors partitioning nodes:
    I1 = []
    for i, node in enumerate(nodes):
        if node.split[long_dim] < centroid[long_dim]:
            I1.append(i)
    I2 = [i for i in range(len(nodes)) if i not in I1]
    assert len(I1) > 0 and len(I2) > 0

    # Swap the index vectors if len(I1) > len(I2):
    n1 = sum(len(nodes[i].get_points()) for i in I1)
    n2 = sum(len(nodes[i].get_points()) for i in I2)
    if n2 < n1:
        I1, I2 = I2, I1

    return I1, I2

def sample_middle_out_butterfly(linOp, rowNodes, colNodes, k, p=8, q):

    m, n = linOp.shape
    M, N = len(rowNodes), len(colNodes)

    # Starting nodes should all come from the same tree and should all
    # be at the same depth
    assert all(node.tree == nodes[0].tree for node in nodes)
    assert all(node.depth == nodes[0].depth for node in nodes)

    assert m == sum(node.get_num_points() for node in rowNodes)
    assert n == sum(node.get_num_points() for node in colNodes)

    i0min = rowNodes[0].get_first_index()
    j0min = colNodes[0].get_first_index()

    assert rowNodes[-1].get_last_index() - i0min == m
    assert colNodes[-1].get_last_index() - j0min == n

    diagBlocks = [[None for _ in range(N)] for _ in range(M)]

    rowFacStreamers = [bf.FacStreamer() for _ in range(M)]
    colFacStreamers = [bf.FacStreamer() for _ in range(N)]

    for colInd, colNode in enumerate(colNodes):
        j0 = colNode.get_first_index() - j0min
        j1 = colNode.get_last_index() - j0min

        # Algorithm 4.1
        Omega = np.zeros((n, k + p))
        Omega[j0:j1]
        Y = linOp@Omega
        Q = np.linalg.qr(Y)[0][:, :k]
        B = Q.T

        assert False
        # U, S, W = ... compute randomized SVD

        colIndPerm = None

        colFacStreamers[colIndPerm].feed(W)

        for rowInd in range(M):
            assert False # TODO: need the column permutation
            diagBlock[rowInd, colIndPerm] = S

        for rowInd, rowNode in enumerate(rowNodes):
            i0 = rowNode.get_first_index() - i0min
            i1 = rowNode.get_last_index() - i0min

            rowFacStreamers.feed(U[i0:i1])

    assert False # extract factors!

    assert False # return!
    # return bf.MatProduct(...)

class DenseLu:
    def __init__(self, A, tol=1e-13):
        m, n = A.shape
        if m != n:
            raise ValueError("A isn't a square matrix")
        del m

        A_dense = np.array(A.to_mat_dense_complex())

        print(A_dense)

        import colorcet as cc
        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm
        plt.figure()
        plt.imshow(complex_to_hsv(A_dense))
        # plt.imshow(abs(A_dense), norm=LogNorm(), cmap=cc.cm.fire_r)
        plt.colorbar()
        plt.show()

        P, L, U = scipy.linalg.lu(A_dense)
        self.P = P
        self.L = L
        self.U = U

class HierarchicalLu:
    '''Prototype class for butterfly fast direct solver for 2D Helmholtz...'''

    DENSE_LU_THRESH = 64

    def __init__(self, A, nodes, depth=0, parent=None):
        print(f'{depth=}, {len(nodes)=}')

        # TODO: handle this case
        assert len(nodes) > 1

        self.A = A
        self.nodes = nodes

        I1, I2 = get_block_inds_for_split(nodes)

        nodes1 = [nodes[i] for i in I1]
        nodes2 = [nodes[i] for i in I2]

        print(f'- {len(nodes1)=}, {len(nodes2)=}')

        A11 = self.A.get_blocks(I1, I1)

        print(f'- {A11.shape=}')

        if max(A11.shape) <= HierarchicalLu.DENSE_LU_THRESH:
            self.A11_lu = DenseLu(A11)
        else:
            self.A11_lu = HierarchicalLu(A11, nodes1, depth=depth+1, parent=self)

        # TODO: this is using the bad old style... but might as well
        # try it and check that it works! Maybe it's fast for a range
        # of problem sizes?

        self.A12 = self.A.get_blocks(I1, I2)
        self.A21 = self.A.get_blocks(I2, I1)

        A22 = self.A.get_blocks(I2, I2)
        A11_refl = bf.MatProduct(self.A21, self.A11_lu, self.A12)
        A11_refl_BF = sample_middle_out_butterfly(A11_refl, nodes2)
        A11_schur = bf.MatDiff(A22, A11_refl_BF)
        if max(A11_schur.shape) <= HierarchicalLu.DENSE_LU_THRESH:
            self.A11_schur_lu = DenseLu(A11_schur)
        else:
            self.A11_schur_lu = HierarchicalLu(A11_schur, nodes2, depth=depth+1, parent=self)

    @staticmethod
    def from_quadtree(A, quadtree, init_level=2):
        return HierarchicalLu(A, quadtree.get_level_nodes(init_level))

    def solve(self, y):
        # y1, y2 = ...
        assert False # XXX

        z1 = self.A11_lu.solve(y1)
        z2 = y2 - self.A21@z1
        x2 = self.A11_schur_lu.solve(z2)
        x1 = z1 - self.A11_lu.solve(self.A12@x2)

        # x = vcat(x1, x2)
        assert False

        return x

############################################################################
# COMMON STUFF FOR REST OF SCRIPT
#

nodes = quadtree.get_level_nodes(2)
block_inds = get_block_inds_for_split(nodes)

############################################################################
# PLOT DISCRETIZED POINTS AND QUADTREE
#

nodes = quadtree.get_level_nodes(2)

if MAKE_PLOTS:
    plt.figure()
    plt.scatter(*np.array(centers).T, s=1, c='r', zorder=1)
    plt.scatter(*np.array(X).T, s=1, c='k', zorder=1)
    plt.gca().add_collection(PatchCollection([
        Rectangle(nodes[i].bbox.xy, nodes[i].bbox.dy, nodes[i].bbox.dx)
        for i in block_inds[0]
    ], facecolor='lightyellow', edgecolor='lightgray', linewidth=1, zorder=0))
    plt.gca().add_collection(PatchCollection([
        Rectangle(nodes[i].bbox.xy, nodes[i].bbox.dy, nodes[i].bbox.dx)
        for i in block_inds[1]
    ], facecolor='lightblue', edgecolor='lightgray', linewidth=1, zorder=0))
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.gca().set_aspect('equal')
    plt.tight_layout()
    plt.show()

############################################################################
# DO A QUICK SANITY CHECK BEFORE TRYING TO BUILD THE HIERARCHICAL LU
# DECOMPOSITION OF A...
#

I1 = np.concatenate([nodes[i].get_inds() for i in block_inds[0]])
I2 = np.concatenate([nodes[i].get_inds() for i in block_inds[1]])

A_dense = np.array(helm.get_kernel_matrix(X, X, N, N))

if MAKE_PLOTS:
    plt.figure()
    plt.imshow(complex_to_hsv(A_dense[I1, :][:, I1]))
    plt.gca().set_aspect('equal')
    plt.show()

A11 = A_dense[I1, :][:, I1]
A12 = A_dense[I1, :][:, I2]
A21 = A_dense[I2, :][:, I1]
A22 = A_dense[I2, :][:, I2]

if MAKE_PLOTS:
    plt.figure()
    plt.imshow(complex_to_hsv(A21@np.linalg.inv(A11)@A12))
    plt.gca().set_aspect('equal')
    plt.show()

def S(mat, normalized=True):
    S = np.linalg.svd(mat)[1]
    if normalized:
        S /= S[0]
    return S

tmp1 = np.linalg.solve(A11, A12)
tmp2 = A21@np.linalg.solve(A11, A12)
tmp3 = A21@np.linalg.inv(A11)

if MAKE_PLOTS:
    plt.figure()
    plt.axhline(y=1e-15, linewidth=1, c='k', zorder=0)
    plt.axvline(x=20, linewidth=1, c='k', zorder=0)
    plt.axvline(x=k+20, linewidth=1, c='k', zorder=0)
    plt.semilogy(S(tmp1), c='r', zorder=1)
    plt.semilogy(S(tmp2), c='b', zorder=1)
    plt.semilogy(S(tmp3), c='g', zorder=1)
    plt.semilogy(S(tmp1[:tmp1.shape[0]//4, :tmp1.shape[1]//4]), c='r', zorder=1)
    plt.semilogy(S(tmp2[:tmp2.shape[0]//4, :tmp2.shape[1]//4]), c='b', zorder=1)
    plt.semilogy(S(tmp3[:tmp3.shape[0]//4, :tmp3.shape[1]//4]), c='g', zorder=1)
    plt.semilogy(S(tmp1[:tmp1.shape[0]//16, :tmp1.shape[1]//16]), c='r', zorder=1)
    plt.semilogy(S(tmp2[:tmp2.shape[0]//16, :tmp2.shape[1]//16]), c='b', zorder=1)
    plt.semilogy(S(tmp3[:tmp3.shape[0]//16, :tmp3.shape[1]//16]), c='g', zorder=1)
    plt.xlim(0, 200)
    plt.ylim(1e-16, 1)
    plt.show()

############################################################################
# COMPUTE THE HIERARCHICAL LU DECOMPOSITION OF A
#

A_lu = HierarchicalLu.from_quadtree(A, quadtree)
