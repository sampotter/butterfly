#!/usr/bin/env python

import numpy as np
import sys

if len(sys.argv) == 1:
    print(f'usage: {sys.argv[0]} <num_points>')
    exit()

num_points = int(sys.argv[1])

theta = np.linspace(0, 2*np.pi, num_points, endpoint=False)

points = np.array([np.cos(theta), np.sin(theta)]).T
normals = points
weights = (theta[1] - theta[0])*np.ones(num_points)
sources = np.array([[0.2, 0.1]])

targets_r = 5
targets_theta = np.linspace(0, 2*np.pi, 2**6, endpoint=False)
targets = targets_r*np.array([np.cos(targets_theta), np.sin(targets_theta)]).T

points.tofile(f'circle{num_points}_points.bin')
normals.tofile(f'circle{num_points}_normals.bin')
weights.tofile(f'circle{num_points}_weights.bin')
sources.tofile(f'circle{num_points}_sources.bin')
targets.tofile(f'circle{num_points}_targets.bin')
