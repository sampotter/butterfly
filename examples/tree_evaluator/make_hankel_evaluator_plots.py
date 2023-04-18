#!/usr/bin/env python

import colorcet as cc
import itertools as it
import matplotlib.pyplot as plt; plt.ion()
import numpy as np
import subprocess

from matplotlib.colors import LogNorm
from numpy.polynomial import Chebyshev, Polynomial
from scipy.special import j0, jn_zeros
from scipy.optimize import brentq

MAKE_PLOTS = True

R0 = 2.4048255576957724
R1 = 20
NUM_POINTS = 1e6
POINTS_TYPE = 'uniform'
TOL = 1e-13
dmin, dmax = 4, 15
kmin, kmax = 2, 15

zeros = np.array([(m - 1/4)*np.pi for m in range(1, 21)])
zeros_fixed = np.array([brentq(j0, x0 - np.pi/8, x0 + np.pi/8, xtol=np.finfo(np.float64).eps) for x0 in zeros])

breakpoints = np.concatenate([[0, 1.2024127788478864], jn_zeros(0, 20), jn_zeros(1, 20)])
breakpoints = np.sort(breakpoints)

plt.figure()
for r0 in zeros_fixed:
    plt.axvline(x=r0, c='k', zorder=1)
plt.plot(np.linspace(0, breakpoints[-1], 501), (np.linspace(0, breakpoints[-1], 501)/np.pi + 1/4) % 1)
plt.show()

def find_degree(r0, r1):
    d = 40
    c = Chebyshev.interpolate(j0, d, [r0, r1]).coef
    return np.where(abs(c/c[0]) < np.finfo(np.float64).resolution)[0][0] - 1

degrees = []
for r0, r1 in zip(breakpoints[:-1], breakpoints[1:]):
    print(r0, r1)
    degrees.append(find_degree(r0, r1))
degree = int(np.median(degrees))

degree = 13
chebs = []
for r0, r1 in zip(breakpoints[:-1], breakpoints[1:]):
    cheb = Chebyshev.interpolate(j0, degree, [r0, r1])
    chebs.append(cheb)

plt.figure()
for cheb in chebs:
    _ = np.linspace(*cheb.domain)
    plt.plot(_, j0(_) - cheb(_), c='k')
plt.show()

plt.figure()
plt.axhline(y=0, c='gray', linewidth=1, zorder=0)
plt.plot(np.linspace(0, breakpoints.max(), int((10*breakpoints.max())//1)), j0(np.linspace(0, breakpoints.max(), int((10*breakpoints.max())//1))), linewidth=2, zorder=1)
plt.scatter(zeros, j0(zeros), zorder=1)
plt.scatter(zeros_fixed, j0(zeros_fixed), c='k', zorder=2)
plt.scatter(breakpoints, j0(breakpoints), marker='x', c='r', zorder=3)
plt.xlim(0, breakpoints.max())
plt.show()

def get_rate(d, k):
    output = subprocess.check_output([f"./test_hankel_evaluator --r0={R0} --r1={R1} --num_points={int(NUM_POINTS//1)} --points_type={POINTS_TYPE} --degree={d} --valence={k} --tol={TOL}"], shell=True, text=True)
    print(output)
    lines = output.split('\n')
    t_build = float(lines[-4].split(':')[-1].split()[0])
    r_sqrt = float(lines[3].split(':')[-1].split()[0])
    r_std = float(lines[5].split(':')[-1].split()[0])
    r_gsl = float(lines[7].split(':')[-1].split()[0])
    r_clenshaw = float(lines[9].split(':')[-1].split()[0])
    r = float(lines[-3].split(':')[-1].split()[0])
    e_max_abs = float(lines[-2].split()[-1])
    return t_build, r_sqrt, r_std, r_gsl, r_clenshaw, r, e_max_abs

D = np.arange(dmin, dmax + 1)

K = np.arange(kmin, kmax + 1)

T_build = np.empty((D.size, K.size))
R_sqrt = np.empty((D.size, K.size))
R_std = np.empty((D.size, K.size))
R_gsl = np.empty((D.size, K.size))
R_clenshaw = np.empty((D.size, K.size))
R = np.empty((D.size, K.size))
E_max_abs = np.empty((D.size, K.size))

for (i, d), (j, k) in it.product(enumerate(D), enumerate(K)):
    print(d, k)
    t_build, r_sqrt, r_std, r_gsl, r_clenshaw, r, e_max_abs = get_rate(d, k)
    T_build[i, j] = t_build
    R_sqrt[i, j] = r_sqrt
    R_std[i, j] = r_std
    R_gsl[i, j] = r_gsl
    R_clenshaw[i, j] = r_clenshaw
    R[i, j] = r
    E_max_abs[i, j] = e_max_abs

iarg, jarg = np.unravel_index(np.argmax(R), R.shape)
print('best combination:')
print(f'- degree = {D[iarg]}')
print(f'- valence = {K[jarg]}')
print(f'- rate = {R[iarg, jarg]/1e6:.1f} Mp/cs')

if MAKE_PLOTS:
    titles_and_values = [
        (r'$t_{build}$', T_build),
        (r'$R_{sqrt}$', R_sqrt),
        (r'$R_{std}$', R_std),
        (r'$R_{gsl}$', R_gsl),
        (r'$R_{clenshaw}$', R_clenshaw),
        (r'$R$', R),
        (r'Max abs error', E_max_abs),
        ('Speedup (vs C standard library)', R/R_std.mean()),
        ('Speedup (vs GSL)', R/R_gsl.mean()),
        ('Tree eval efficiency (wrt Clenshaw)', R/R_clenshaw),
        ('C std efficiency (wrt Clenshaw)', R_std/R_clenshaw),
        ('GSL efficiency (wrt Clenshaw)', R_gsl/R_clenshaw)
    ]
    for title, values in titles_and_values:
        plt.figure()
        kwargs = {
            'extent': [kmin - 0.5, kmax + 0.5, dmin - 0.5, dmax + 0.5],
            'cmap': cc.cm.bmy,
            'interpolation': 'none'
        }
        if title == 'Max abs error':
            kwargs['norm'] = LogNorm(vmin=1e-16, vmax=1e-10)
        plt.imshow(np.flipud(values), **kwargs)
        plt.colorbar()
        plt.scatter(K[jarg], D[iarg], marker='x', s=75, c='k')
        plt.gca().set_xticks(K)
        plt.gca().set_yticks(D[::-1])
        plt.xlabel('$k$')
        plt.ylabel('$d$')
        plt.title(title)
        plt.gca().set_aspect('equal')
        plt.tight_layout()
        plt.show()
