#!/usr/bin/env python

import numpy as np
import sys

if len(sys.argv) == 1:
    print(f'usage: {sys.argv[0]} <num_points>')
    exit()

num_points = int(sys.argv[1])
num_targets = 5

# set up ellipse

x0, y0 = 1, 0.2 # point source coordinates
a, b = 1.5, 2 # radii of ellipse
r_eval = 2*max(a, b) # evaluation radius

T = np.linspace(0, 2*np.pi, num_points, endpoint=False)
h = 2*np.pi/num_points

x = a*np.cos(T) + x0
y = b*np.sin(T) + y0
Q = np.array([x, y]).T

dx, d2x = -a*np.sin(T), -a*np.cos(T)
dy, d2y =  b*np.cos(T), -b*np.sin(T)
dsdt = np.sqrt(dx**2 + dy**2)
tx, ty = dx/dsdt, dy/dsdt
nx, ny = ty, -tx # outward-facing

kappa = (dx*d2y - dy*d2x)/dsdt**3

# write data to binary files

points = np.array([x, y]).T
normals = np.array([nx, ny]).T
weights = h*dsdt
sources = np.array([[x0, y0]])

targets_theta = np.linspace(0, 2*np.pi, num_targets, endpoint=False)
targets = r_eval*np.array([np.cos(targets_theta), np.sin(targets_theta)]).T

points.tofile(f'ellipse{num_points}_points.bin')
normals.tofile(f'ellipse{num_points}_normals.bin')
weights.tofile(f'ellipse{num_points}_weights.bin')
sources.tofile(f'ellipse{num_points}_sources.bin')
targets.tofile(f'ellipse{num_points}_targets.bin')
