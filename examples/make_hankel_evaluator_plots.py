import colorcet as cc
import itertools as it
import matplotlib.pyplot as plt; plt.ion()
import numpy as np
import subprocess

from matplotlib.colors import LogNorm

R0 = 0
R1 = 500
NUM_POINTS = 1e6
POINTS_TYPE = 'uniform'
TOL = 1e-13

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

dmin, dmax = 4, 10
D = np.arange(dmin, dmax + 1)

kmin, kmax = 2, 10
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
