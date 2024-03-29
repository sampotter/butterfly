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
import itertools as it
import numpy as np
import scipy.linalg

randn = np.random.randn

# Shortcut for Hermitian tranpose:
H = lambda _: _.conj().T

import butterfly as bf

from util import complex_to_hsv

if MAKE_PLOTS:
    import colorcet as cc
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
h = 0.005 # point spacing
KR_order = 6 # Kapur-Rokhlin quadrature order
eps = 1e-13 # Relative error tolerance for hierarchial LU decomposition

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

if MAKE_PLOTS:
    plt.figure()
    # quadtree.plot_node_boxes()
    for node in quadtree.nodes:
        bbox = node.bbox
        x0, x1, y0, y1 = bbox.xmin, bbox.xmax, bbox.ymin, bbox.ymax
        plt.plot([x0, x1, x1, x0, x0], [y0, y0, y1, y1, y0], c='r', zorder=2)
    plt.scatter(*np.array(X).T, s=2.5, c='k', zorder=1)
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.gca().set_aspect('equal')
    plt.show()

# get permuted copy of quadrature weights
rev_perm = quadtree.perm.get_reverse()
W_perm = W.copy()
W_perm.permute(rev_perm)

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

def hlu_recursive_base_case(A, nodes, p):
    return A.shape[0] <= HierarchicalLu.DENSE_LU_THRESH \
        or len(nodes) == 1 \
        or any(_.get_num_points() <= p for _ in nodes)

def minimize_rank_est(X1, X2, k, eps):
    r1min = lambda x1, y1: np.sqrt(np.sum(np.subtract(X1, (x1, y1))**2, axis=1)).max()
    r2min = lambda x2, y2: np.sqrt(np.sum(np.subtract(X2, (x2, y2))**2, axis=1)).max()
    R = lambda x1, y1, x2, y2: np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
    f = lambda x1, y1, x2, y2, r1, r2: r1*r2/(R(x1, y1, x2, y2) - r1 - r2)
    g1 = lambda x1, y1, r1: r1 - r1min(x1, y1)
    g2 = lambda x2, y2, r2: r2 - r2min(x2, y2)
    g3 = lambda x1, y1, x2, y2, r1, r2: R(x1, y1, x2, y2) - r1 - r2
    p1, p2 = X1.mean(0), X2.mean(0)
    _ = scipy.optimize.minimize(
        lambda _: f(*_),
        (*p1, *p2, r1min(*p1), r2min(*p2)),
        bounds=[
            (None, None), # x1
            (None, None), # y1
            (None, None), # x2
            (None, None), # y2
            (0, None),    # r1
            (0, None)],   # r2
        constraints=[
            {'type': 'ineq', 'fun': lambda _: g1(_[0], _[1], _[4])},
            {'type': 'ineq', 'fun': lambda _: g2(_[2], _[3], _[5])},
            {'type': 'ineq', 'fun': lambda _: g3(*_)}],
        options={'maxiter': 1000, 'ftol': 1e-3})
    if _.success:
        p1opt, p2opt, r1opt, r2opt = _.x[:2], _.x[2:4], _.x[4], _.x[5]

        # fValue = f(*p1opt, *p2opt, r1min(*p1opt), r2min(*p2opt))
        assert _.fun > 0
        rank = min(X1.shape[0], X2.shape[0], k*_.fun + np.log10(1/eps))
        print(f'{rank = }')

        # import matplotlib.pyplot as plt
        # plt.figure()
        # plt.scatter(*X1.T, s=2, c='r', zorder=1)
        # plt.scatter(*X2.T, s=2, c='b', zorder=1)
        # plt.scatter(*p1opt, s=10, c='r', zorder=2)
        # plt.scatter(*p2opt, s=10, c='b', zorder=2)
        # plt.scatter([p1opt[0], p2opt[0]], [p1opt[1], p2opt[1]], c='k', zorder=0)
        # plt.gca().add_patch(plt.Circle(p1opt, r1opt, edgecolor='r', facecolor='none'))
        # plt.gca().add_patch(plt.Circle(p2opt, r2opt, edgecolor='b', facecolor='none'))
        # plt.gca().set_aspect('equal')
        # plt.show()

        return p1opt, p2opt, r1opt, r2opt, rank

# def rank_for_node_is_small_enough(node, X2, k, eps, p):
#     if node.get_num_points() <= p:
#         return True

#     X1 = np.array(node.get_points())

#     _ = minimize_rank_est(X1, X2, k, eps)
#     if _ is None:
#         return False
#     pX, pY, rX, rY, popt = _
#     print(f'{popt = :.2f}')

#     # if MAKE_PLOTS:
#     import matplotlib.pyplot as plt
#     plt.figure()
#     plt.scatter(*X1.T, s=5, c='b')
#     plt.scatter(*X2.T, s=5, c='r')
#     plt.scatter(*pX, s=10, c='b', marker='x')
#     plt.scatter(*pY, s=10, c='r', marker='x')
#     plt.gca().add_patch(plt.Circle(pX, rX, edgecolor='k', facecolor='none'))
#     plt.gca().add_patch(plt.Circle(pY, rY, edgecolor='k', facecolor='none'))
#     plt.gca().set_aspect('equal')
#     plt.show()

#     return popt <= p

def estimate_rank(srcNode, tgtNodes, k, eps, p):
    numSrcPoints = srcNode.get_num_points()
    numTgtPoints = sum(_.get_num_points() for _ in tgtNodes)
    maxNumPoints = max(numSrcPoints, numTgtPoints)
    if maxNumPoints <= p:
        return maxNumPoints

    Xsrc = np.array(srcNode.get_points())

    def get_results(tgtNode):
        return minimize_rank_est(Xsrc, np.array(tgtNode.get_points()), k, eps)

    def get_p(res): return res[-1]

    cutNodes = tgtNodes.copy()
    results = [get_results(_) for _ in cutNodes]

    i = 0
    while i < len(cutNodes):
        cutNode, res = cutNodes[i], results[i]
        childResults = [get_results(_) for _ in cutNode.children]

        # If the optimization failed for the current node, we push
        # down to its children and keep going:
        if res is None:
            cutNodes = cutNodes[:i] + cutNode.children + cutNodes[(i + 1):]
            results = results[:i] + childResults + results[(i + 1):]
            continue

        # If we were able to optimize for the current node, we should
        # be able to optimize for all of its children:
        try:
            assert all(_ is not None for _ in childResults)
        except:
            for j, _ in enumerate(cutNode.children):
                if childResults[j] is None:
                    import ipdb; ipdb.set_trace()
                    get_results(_)

        # If the child nodes together have a lower rank than the
        # current node, splice them in and move on:
        if sum(get_p(_) for _ in childResults) <= get_p(res):
            cutNodes = cutNodes[:i] + cutNode.children + cutNodes[(i + 1):]
            results = results[:i] + childResults + results[(i + 1):]
            continue

        # The current node minimizes the rank, so we can keep it
        i += 1

    # import matplotlib.pyplot as plt
    # plt.figure()
    # plt.scatter(*Xsrc.T, s=2, c='k', zorder=1)
    # c = ['r', 'g', 'orange', 'yellow', 'magenta', 'brown']
    # j = 0
    # for cutNode, res in zip(cutNodes, results):
    #     p1opt, p2opt, r1opt, r2opt = res[:4]
    #     plt.scatter(*np.array(cutNode.get_points()).T, s=2, c=c[j], zorder=1)
    #     plt.scatter(*p1opt, s=10, c='k', zorder=2)
    #     plt.scatter(*p2opt, s=10, c=c[j], zorder=2)
    #     plt.plot([p1opt[0], p2opt[0]], [p1opt[1], p2opt[1]], c='k', zorder=0)
    #     plt.gca().add_patch(plt.Circle(p1opt, r1opt, edgecolor='k', facecolor='none'))
    #     plt.gca().add_patch(plt.Circle(p2opt, r2opt, edgecolor=c[j], facecolor='none'))
    #     j += 1
    # plt.gca().set_aspect('equal')
    # plt.tight_layout()
    # plt.show()

    return sum(get_p(_) for _ in results)

def rank_for_all_nodes_is_small_enough(nodes, X2, k, eps, p):
    return all(rank_for_node_is_small_enough(_, X2, k, eps, p) for _ in nodes)

def get_next_tree_level(nodes):
    assert len(nodes) == 1 or all(nodes[0].depth == _.depth for _ in nodes[1:])
    childNodes = []
    for node in nodes:
        childNodes.extend(node.children)
    return childNodes

def get_common_parent(nodeSpan):
    assert all(_.depth == nodeSpan[0].depth for _ in nodeSpan[1:])
    parents = {_.parent for _ in nodeSpan}
    while len(parents) > 1:
        parents = {_.parent for _ in parents}
    return list(parents)[0]

def estimate_compression_rate(reflNodes, tgtNodes, k, eps, p):
    sizes = [_.get_num_points() for _ in reflNodes]
    rankEsts = [estimate_rank(_, tgtNodes, k, eps, p) for _ in reflNodes]
    blockRanks = [[min(r1, r2) for r2 in rankEsts] for r1 in rankEsts]
    compressedSize = sum(
        sum((n1 + n2)*min(r1, r2) for n1, r1 in zip(sizes, rankEsts))
        for n2, r2 in zip(sizes, rankEsts)
    )
    uncompressedSize = sum(sizes)**2
    return uncompressedSize/compressedSize

def get_refl_nodes(srcNodes, tgtNodes, k, eps, p):
    # # Compute the minimum number of nodes we require in a level
    # sqrtN = int(np.ceil(np.sqrt(srcParent.get_num_points())))
    # Xtgt = np.concatenate([_.get_points() for _ in tgtNodes], axis=0)

    def get_R(_):
        return estimate_compression_rate(_, tgtNodes, k, eps, p)

    reflNodes = srcNodes.copy()
    nextReflNodes = get_next_tree_level(reflNodes)
    R, nextR = get_R(reflNodes), get_R(nextReflNodes)
    print(f'{R = }')
    print(f'R = {nextR}')
    while nextR > R:
        reflNodes = nextReflNodes
        nextReflNodes = get_next_tree_level(reflNodes)
        R, nextR = nextR, get_R(nextReflNodes),
        print(f'R = {nextR}')

    return reflNodes

def zrandn(m, n):
    return randn(m, n)/np.sqrt(2) + 1j*randn(m, n)

def sample_middle_out_butterfly(linOp, rowNodes, colNodes, eps, p, q):
    '''Sample a middle-out butterfly factorization of `linOp`.

    The middle blocks are initially compressed using a randomized SVD.

    Args:
        linOp: the linear operator compress (must support linOp.shape and linOp@x)
        rowNodes: the tree nodes giving the row indices of the starting blocks in linOp
        colNodes: like rowNodes, but for column blocks
        eps: the target tolerance for the low-rank approximation
        p: the target rank of the low-rank approximation
        q: the oversampling parameter for the randomized SVD

    Returns:
        bf.MatProduct: the butterfly factorization

    '''

    m, n = linOp.shape
    M, N = len(rowNodes), len(colNodes)

    # Things get a little crazy with the starting low-rank
    # factorizations below if we allow the row or column nodes to
    # contain fewer than `p` points... So, just ensure that there are
    # at least `p` points here. We could look into relaxing this later
    # but probably not worth the complication.
    # assert all(node.get_num_points() >= p for node in rowNodes)
    # assert all(node.get_num_points() >= p for node in colNodes)
    if not all(node.get_num_points() >= p for node in rowNodes) or \
       not  all(node.get_num_points() >= p for node in colNodes):
        print('ack!')
        B = np.array(linOp.to_mat_dense_complex())
        import ipdb; ipdb.set_trace()

    # Starting nodes should all come from the same tree and should all
    # be at the same depth
    for nodes in [rowNodes, colNodes]:
        assert all(node.tree == nodes[0].tree for node in nodes)
        assert all(node.depth == nodes[0].depth for node in nodes)

    assert m == sum(node.get_num_points() for node in rowNodes)
    assert n == sum(node.get_num_points() for node in colNodes)

    rowSubtree = bf.Tree.from_nodes(rowNodes)
    colSubtree = bf.Tree.from_nodes(colNodes)

    # Now I need to set up two trees with the same structure as
    # `rowSubtree` and `colSubtree`, but which are just plain "index
    # trees" and whose leaves all contain `p` nodes.
    #
    #(Again, see the comment above about trying to make the starting
    # low-rank factorizations smaller than `p`---would be a lot of
    # added complication, probably for little benefit. Can consider
    # this idea later as an optimization...)

    rowIndexSubtree = bf.Tree.for_middle_fac(rowSubtree, p)
    colIndexSubtree = bf.Tree.for_middle_fac(colSubtree, p)

    if MAKE_PLOTS:
        plt.figure()
        ax = plt.gca()
        bf.Quadtree.from_tree(rowSubtree).plot_node_boxes(ax)
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.set_aspect('equal')
        plt.tight_layout()
        plt.show()

    i0min = rowNodes[0].i0
    j0min = colNodes[0].i0

    kwargs = dict(tol=eps, minNumRows=p, minNumCols=p)

    rowFacStreamers, rowPerms = [], []
    for a in range(M):
        subtree, perm = bf.Tree.from_node(rowNodes[a])
        facStreamer = bf.FacStreamer.from_trees(subtree, colIndexSubtree, **kwargs)
        rowFacStreamers.append(facStreamer)
        rowPerms.append(perm)

    colFacStreamers, colPerms = [], []
    for b in range(N):
        subtree, perm = bf.Tree.from_node(colNodes[b])
        facStreamer = bf.FacStreamer.from_trees(subtree, rowIndexSubtree, **kwargs)
        colFacStreamers.append(facStreamer)
        colPerms.append(perm)

    # Set up the nonzero part of the left and right random test
    # matrices. We need to precompute these so that we can reuse them
    # later when we solve the least-squares problems for the middle
    # blocks (the B's).
    #
    # TODO: might not actually need to store Omega...
    OmegaBlocks = [zrandn(_.get_num_points(), p + q) for _ in colNodes]
    OmegaTildeBlocks = [zrandn(_.get_num_points(), p + q) for _ in rowNodes]

    ABlocks = np.empty((M, N), dtype=object)
    BBlocks = np.empty((M, N), dtype=object)

    # Sample the range of each block column, stream the left butterfly
    # factors we pick up, and collect the system matrix for each least
    # square problem we need to solve to get the middle factors.
    for b, colNode in enumerate(colNodes):
        j0 = colNode.i0 - j0min
        j1 = colNode.i1 - j0min

        Omega = np.zeros((n, p + q), np.complex128)
        Omega[j0:j1] = OmegaBlocks[b]
        Y = np.array(linOp@Omega)

        for a, rowNode in enumerate(rowNodes):
            i0 = rowNode.i0 - i0min
            i1 = rowNode.i1 - i0min

            Q, S = np.linalg.svd(Y[i0:i1], full_matrices=False)[:2]
            assert S.size <= p or S[p] <= eps*S[0]
            Q = np.ascontiguousarray(Q[:, :p])
            rowFacStreamers[a].feed(Q)

            ABlocks[a, b] = H(OmegaTildeBlocks[a])@Q

    assert all(_.is_done() for _ in rowFacStreamers)

    # Sample the range of each block row, stream the right butterfly
    # factors, and collect the load matrices for each of the least
    # squares problems.
    for a, rowNode in enumerate(rowNodes):
        i0 = rowNode.i0 - i0min
        i1 = rowNode.i1 - i0min

        # YTilde = linOp^H * OmegaTilde -- but this is a very
        # inefficient way to do this...
        OmegaTilde = np.zeros((m, p + q), np.complex128)
        OmegaTilde[i0:i1] = OmegaTildeBlocks[a]
        OmegaTildeH = np.ascontiguousarray(H(OmegaTilde))
        YTilde = H(np.array(linOp.__rmatmul__(OmegaTildeH)))

        for b, colNode in enumerate(colNodes):
            j0 = colNode.i0 - j0min
            j1 = colNode.i1 - j0min

            QTilde, STilde = np.linalg.svd(YTilde[j0:j1], full_matrices=False)[:2]
            assert STilde.size <= p or STilde[p] <= eps*STilde[0]
            QTilde = np.ascontiguousarray(QTilde[:, :p])
            colFacStreamers[b].feed(QTilde)

            BBlocks[a, b] = H(YTilde[j0:j1])@QTilde

    assert all(_.is_done() for _ in colFacStreamers)

    # Solve each of the least squares problems for the middle factors:
    MiddleBlocks = np.empty((M, N), dtype=object)
    for a, b in it.product(range(M), range(N)):
        _ = np.linalg.lstsq(ABlocks[a, b], BBlocks[a, b], rcond=None)
        l2_error = np.sqrt(_[1][-1])
        if l2_error > eps:
            print(f'warning: error in middle block ({l2_error:.1e}) exceeds tolerance ({eps:.1e})')
        MiddleBlocks[a, b] = _[0]

    # Compute index offsets for the middle blocks:
    I0 = np.empty((M, N), dtype=int); I0[...] = -1
    J0 = np.empty((M, N), dtype=int); J0[...] = -1
    i0 = 0
    for a in range(M):
        j0 = sum(_.shape[1] for _ in MiddleBlocks[0, :a])
        for b in range(N):
            I0[a, b] = i0
            J0[a, b] = j0
            i0 += MiddleBlocks[a, b].shape[0]
            j0 += sum(_.shape[1] for _ in MiddleBlocks[b, :])

    # Get the shape of the middle factor:
    middleFactorShape = (
        sum(_.shape[0] for _ in  MiddleBlocks.ravel()),
        sum(_.shape[1] for _ in  MiddleBlocks.ravel())
    )

    # Build a list of the middle blocks with their indices (each entry
    # is a tuple of the form (i0, j0, block)):
    indexedMiddleBlocks = list(zip(I0.ravel(), J0.ravel(), MiddleBlocks.ravel()))

    if MAKE_PLOTS:
        plt.figure()
        rects = []
        for i0, j0, middleBlock in indexedMiddleBlocks:
            m, n = middleBlock.shape
            rect = Rectangle((i0, j0), m, n)
            rects.append(rect)
        plt.gca().add_collection(PatchCollection(rects, facecolor='cyan', edgecolor='k', linewidth=1))
        plt.xlim(0, middleFactorShape[0])
        plt.ylim(0, middleFactorShape[1])
        plt.xlabel('$i$')
        plt.ylabel('$j$')
        plt.gca().set_aspect('equal')
        plt.show()

    # Set up factors:
    leftFactor = bf.MatBlockDiag.from_blocks([_.toMat() for _ in rowFacStreamers], bf.Policy.View)
    middleFactor = bf.MatBlockCoo.from_indexed_blocks(middleFactorShape, indexedMiddleBlocks, bf.Policy.View)
    rightFactor = bf.MatBlockDiag.from_blocks([_.toMat() for _ in colFacStreamers], bf.Policy.View)
    rightFactor.transpose()

    # Build and return factorization:
    return bf.MatProduct.from_factors(leftFactor, middleFactor, rightFactor)

class DenseLu(bf.MatPython):
    def __init__(self, A, tol=1e-13):
        m, n = A.shape
        if m != n:
            raise ValueError("A isn't a square matrix")

        A_dense = np.array(A.to_mat_dense_complex())
        if np.linalg.matrix_rank(A_dense) < n:
            import ipdb; ipdb.set_trace()
            raise ValueError("A is singular")

        super().__init__(n, n)

        P, L, U = scipy.linalg.lu(A_dense)
        self.P = P
        self.L = L
        self.U = U

    def _Mul(self, X):
        _ = self.P.T@X
        _ = scipy.linalg.solve_triangular(self.L, _, lower=True, unit_diagonal=True)
        _ = scipy.linalg.solve_triangular(self.U, _)
        return _

    def _Rmul(self, X):
        _ = scipy.linalg.solve_triangular(self.U.T, X.T, lower=True)
        _ = scipy.linalg.solve_triangular(self.L.T, _, unit_diagonal=True)
        _ = self.P@_
        return _.T

class HierarchicalLu(bf.MatPython):
    '''Prototype class for butterfly fast direct solver for 2D Helmholtz...'''

    DENSE_LU_THRESH = 512

    def __init__(self, A, nodes, k, eps, p=20, q=4, depth=0, parent=None):
        m, n = A.shape
        if m != n:
            raise ValueError("A isn't a square matrix")

        super().__init__(n, n)

        self.A = A
        self.nodes = nodes

        I1, I2 = get_block_inds_for_split(nodes)
        nodes1 = [nodes[i] for i in I1]
        nodes2 = [nodes[i] for i in I2]
        assert nodes1 and nodes2

        # Get diagonal blocks:
        A11 = self.A.get_blocks(I1, I1)
        A22 = self.A.get_blocks(I2, I2)

        # We need to hang on to A11 for the solve:
        self.A11 = A11

        if hlu_recursive_base_case(A11, nodes1, p):
            if len(nodes1) == 1:
                print(f'''warning: found leaf node with more points (== {nodes1[0].get_num_points()}) than HierarchicalLu.DENSE_LU_THRESH (== {HierarchicalLu.DENSE_LU_THRESH}) while recursively factorizing leading subblock''')
            print('  ' * (depth + 1), end='')
            print(f'DenseLu(n={A11.shape[0]}) [leading block]')
            self.A11_lu = DenseLu(A11)
        else:
            print('  ' * (depth + 1), end='')
            print(f'HierarchicalLu(n={A11.shape[0]}) [leading block]')
            self.A11_lu = HierarchicalLu(A11, nodes1, k, eps, depth=depth+1, parent=self)

        # TODO: this is using the bad old style... but might as well
        # try it and check that it works! Maybe it's fast for a range
        # of problem sizes?
        #
        # (specifically, it might be possible to speed this up
        # asymptotically if we butterfly the LU factors (or, rather,
        # butterfly A21 A11^-1 and A11^-1 A12... empirically it seems
        # like these can be butterflied...)

        self.A12 = self.A.get_blocks(I1, I2)
        self.A21 = self.A.get_blocks(I2, I1)

        # Set up reflector:
        A11_refl = bf.MatProduct.from_factors(self.A21, self.A11_lu, self.A12)

        if MAKE_PLOTS:
            plt.figure()
            plt.imshow(np.real(np.array(A11_refl.to_mat_dense_complex())), cmap=cc.cm.gouldian)
            plt.colorbar()
            plt.show()

        A11_refl_nodes = get_refl_nodes(nodes2, nodes1, k, eps, p)
        if A11_refl_nodes:
            A11_refl_BF = sample_middle_out_butterfly(
                A11_refl, A11_refl_nodes, A11_refl_nodes, eps, p, q)
            A11_schur = bf.MatDiff(A22, A11_refl_BF)
        else:
            if A11_refl.shape[0] > HierarchicalLu.DENSE_LU_THRESH:
                print(f'warning: failed to BF reflector for node with {n} points')
                import ipdb; ipdb.set_trace()
                get_refl_nodes(nodes2, nodes1, k, eps, p)
            A11_schur = bf.MatDiff(A22, A11_refl)

        # Continue factorization recursively:
        if hlu_recursive_base_case(A11_schur, nodes2, p):
            if len(nodes2) == 1:
                print(f'''warning: found leaf node with more points (== {nodes2[0].get_num_points()}) than HierarchicalLu.DENSE_LU_THRESH (== {HierarchicalLu.DENSE_LU_THRESH}) while recursively factorizing Schur complement''')
            print('  ' * (depth + 1), end='')
            print(f'DenseLu(n={A11_schur.shape[0]}) [Schur complement]')
            self.A11_schur_lu = DenseLu(A11_schur)
        else:
            print('  ' * (depth + 1), end='')
            print(f'HierarchicalLu(n={A11_schur.shape[0]}) [Schur complement]')
            self.A11_schur_lu = \
                HierarchicalLu(A11_schur, nodes2, k, eps, depth=depth+1, parent=self)

    @staticmethod
    def from_quadtree(A, quadtree, k, eps, init_level=2):
        print(f'HierarchicalLu.from_quadtree(n={A.shape[0]})')
        return HierarchicalLu(A, quadtree.get_level_nodes(init_level), k, eps)

    @property
    def i0(self):
        return 0

    @property
    def i1(self):
        return self.A11_lu.shape[0]

    @property
    def i2(self):
        return self.shape[0]

    @property
    def j0(self):
        return 0

    @property
    def j1(self):
        return self.A11_lu.shape[1]

    @property
    def j2(self):
        return self.shape[1]

    def _Mul(self, X):
        assert self.shape[1] == X.shape[0]
        X1, X2 = X[self.j0:self.j1], X[self.j1:self.j2]
        Z1 = self.A11_lu@X1
        Z2 = X2 - self.A21@Z1
        Y2 = self.A11_schur_lu@Z2
        Y1 = self.A11@X1 - self.A11_lu@(self.A12@Y2)
        Y = np.concatenate([Y1, Y2], axis=0)
        assert Y.shape[0] == self.shape[0]
        assert Y.shape[1] == X.shape[1]
        return Y

    def _Rmul(self, X):
        print('HierarchicalLu._Rmul')
        assert False

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

print('assembling hierarchical LU decomposition:')
A_lu = HierarchicalLu.from_quadtree(A, quadtree, k, eps)
