#!/usr/bin/env python

import itertools as it
import numpy as np
import subprocess
import time

from pathlib import Path

MIN_POINTS = 10

# The K grid for tests. We set up several uniforms grids at different
# scales. The first one from k = 1 to 9, the second from k = 10 to 90,
# and the final one from k = 100 to 1000 (with the obvious spacing).
K = np.concatenate([
    np.linspace(1.0, 9.0, 9),
    np.linspace(10.0, 90.0, 9),
    np.linspace(100.0, 1000.0, 10)
])

# The H grid for tests. We set it up to match K so that each value
# will give roughly 10 points per wavelength for that wavenumber.
H = 2*np.pi/(MIN_POINTS*K)

R = np.array([0.025, 0.05,  0.1,  0.2,  0.3,  0.4, 0.5])[::-1]
A = np.array([0.005, 0.01, 0.02, 0.04, 0.06, 0.08, 0.1])[::-1]
B = np.array([ 0.01, 0.02, 0.04, 0.08, 0.12, 0.16, 0.2])[::-1]

def get_lambda(k):
    return 2*np.pi/k

def get_points_per_wavelength(k, h):
    return get_lambda(k)/h

def get_points_per_ellipse(a, b, h):
    return (a + b)/h

if __name__ == '__main__':
    results_dir_path = Path('results')
    results_dir_path.mkdir(exist_ok=True)

    results_raw_dir_path = results_dir_path/'raw'
    results_raw_dir_path.mkdir(exist_ok=True)

    for k, (r, a, b), h in it.product(K, zip(R, A, B), H):
        ppw = get_points_per_wavelength(k, h)
        ppe = get_points_per_ellipse(a, b, h)
        if ppw < MIN_POINTS or ppe < MIN_POINTS:
            continue

        output_path = results_raw_dir_path/f'k{k}_r{r}_a{a}_b{b}_h{h:0.1E}.txt'
        if output_path.exists():
            continue

        print(f'{k=}, {r=}, {a=}, {b=}, {h=:0.1E}, {ppw=:0.1f}, {ppe=:0.1f}...', end='', flush=True)

        t0 = time.time()

        cmd = ['./multiple_scattering', f'-k{k}', f'-r{r}', f'-a{a}', f'-b{b}', f'-h{h}',
               '--doButterflyPreconditionedGmres', '--doFmmPreconditionedGmres']

        proc = subprocess.run(cmd, capture_output=True)

        tf = time.time()

        print(f' done [{tf - t0:.02f}s]')

        with open(output_path, 'w') as f:
            print(proc.stdout.decode('ascii'), file=f)
