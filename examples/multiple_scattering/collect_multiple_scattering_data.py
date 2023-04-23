#!/usr/bin/env python

import itertools as it
import numpy as np
import subprocess
import time

from pathlib import Path

MIN_POINTS = 10
MAX_POINTS = 250000

K = np.logspace(0, 3, 13)

# The H grid for tests. We set it up to match K so that each value
# will give roughly 10 points per wavelength for that wavenumber.
# H = 2*np.pi/(MIN_POINTS*K)

def get_H():
    k = 0
    while True:
        h = (2*np.pi/MIN_POINTS)*10**(-k/4)
        yield h
        k += 1

R = np.array([ 0.1,  0.2,  0.3,  0.4, 0.5])[::-1]
A = R/5
B = 2*A

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

    for k, (r, a, b) in it.product(K, zip(R, A, B)):
        for h in get_H():
            ppw = get_points_per_wavelength(k, h)
            ppe = get_points_per_ellipse(a, b, h)
            if ppw < MIN_POINTS or ppe < MIN_POINTS:
                continue

            output_path = results_raw_dir_path/f'k{k}_r{r}_a{a}_b{b}_h{h:0.1E}.txt'
            if output_path.exists():
                continue

            # run multiple_scattering just to set up the
            # discretization of the domain to collect some stats
            cmd = ['./multiple_scattering', f'-k{k}', f'-r{r}', f'-a{a}', f'-b{b}', f'-h{h}']
            proc = subprocess.run(cmd, capture_output=True)
            lines = proc.stdout.decode('ascii').split('\n')

            # bail if there are too many points in the domain (will
            # run into memory issues)
            line = next(_ for _ in lines if 'Setting up discretization...' in _)
            num_points = int(line.split()[-3])
            if num_points > MAX_POINTS:
                break

            # bail if there are too many points per ellipse (block
            # Jacobi preconditioner will be way too expensive)
            line = next(_ for _ in lines if 'Setting up problem geometry...' in _)
            num_ellipses = int(line.split()[-3])
            if num_ellipses < np.sqrt(num_points)/2:
                break
            #
            print(f'{k=}, {r=}, {a=}, {b=}, {h=:0.1E}:')
            print(f'- {ppw=:0.1f}')
            print(f'- {ppe=:0.1f}')
            print(f'- {num_points=}')

            t0 = time.time()

            cmd = ['./multiple_scattering', f'-k{k}', f'-r{r}', f'-a{a}', f'-b{b}', f'-h{h}',
                   '--doButterflyPreconditionedGmres', '--doFmmPreconditionedGmres']

            proc = subprocess.run(cmd, capture_output=True)

            tf = time.time()

            print(f'- finished: elapsed time = [{tf - t0:.02f}s]')

            with open(output_path, 'w') as f:
                print(proc.stdout.decode('ascii'), file=f)
