#!/usr/bin/env python3

import sys; sys.path.insert(-1, '../../wrappers/python')

import argparse
import matplotlib.pyplot as plt
import numpy as np

import butterfly as bf

# simulation parameters
xmin, xmax = -1, 1
ymin, ymax = -1, 1
xsrc = bf.Points2.from_point((-0.5, 2))
k = 50
r = 0.4 # spacing between ellipses
a, b = 0.04, 0.08 # range of ellipse semi axes
h = 0.02 # point spacing
KR_order = 6 # Kapur-Rokhlin quadrature order

# seed bf's PRNG
bf.seed(0)

# use Poisson disk sampling to generate ellipse centers
bbox = bf.Bbox2(xmin, xmax, ymin, ymax)
centers = bf.Points2.sample_poisson_disk(bbox, r)
num_ellipses = len(centers)

# randomly sample ellipses
ellipses = []
for i in range(num_ellipses):
    a, b = np.random.uniform((a, b), size=2)
    theta = np.random.uniform(2*np.pi)
    ellipses.append(bf.Ellipse(max(a, b), min(a, b), centers[i], theta))

# discretize the ellipses to get the problem geometry
X, W, N = bf.Points2(), bf.RealArray(), bf.Vectors2()
for ellipse in ellipses:
    n = int(np.ceil(ellipse.perimeter/h))
    X_, _, N_, W_ = ellipse.sample_linspaced(n)
    X.extend(X_)
    N.extend(N_)
    W.extend(W_)
print(f'set up discretization with {len(X)} points')

# build quadtree on points w/ normals
quadtree = bf.Quadtree(X, N)
rev_perm = quadtree.perm.get_reverse()

# compute the incident field
phi_in = bf.Helm2.get_kernel_matrix(xsrc, k, bf.LayerPot.Sp, Xtgt=X, Ntgt=N)

# set up multilevel butterfly factorization
A_BF = bf.FacHelm2.make_multilevel(quadtree, k, bf.LayerPot.Sp)
A_BF.apply_KR_correction(KR_order)
A_BF.scale_cols(w[rev_perm])
A_BF += bf.MatIdentity()/2
