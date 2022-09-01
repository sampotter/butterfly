#!/usr/bin/env python

import numpy as np
import sys

if len(sys.argv) == 1:
    print(f'usage: {sys.argv[0]} <num_points>')
    exit()

num_points = int(sys.argv[1])

theta = np.linspace(0, 2*np.pi, num_points, endpoint=False)

points = np.array([np.cos(theta), np.sin(theta)]).T

points.tofile('circle_with_%d_points.bin' % theta.size)
