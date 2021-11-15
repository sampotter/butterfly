import numpy as np

theta = np.linspace(0, 2*np.pi, 2**14, endpoint=False)

points = np.array([np.cos(theta), np.sin(theta)]).T

points.tofile('circle_with_%d_points.bin' % theta.size)
