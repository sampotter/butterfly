import matplotlib.pyplot; plt.ion()
import numpy as np

# quick helper script to remind me how to sample points on curves with
# a given density---these examples show "uniform" density
# (i.e. roughly sampling from an arc length parametrized curve) and
# "1/kappa" (i.e. density matches the curvature)

a, b = 3, 0.5
p = np.pi*(a + b)*(1 + ((a - b)/(a + b))**2/4) # estimate perimeter

# arc-length sampled points on an ellipse

h = 0.1
N = int(p//h)
if N % 2 == 1: N += 1
T = np.linspace(0, 2*np.pi, N + 1)
x, y = a*np.cos(T), b*np.sin(T)
dx = x[1:] - x[:-1]
dy = y[1:] - y[:-1]
D = np.concatenate([[0], np.cumsum(np.sqrt(dx**2 + dy**2))])
T1 = np.interp(np.linspace(D[0], D[-1], T.size), D, T)
x1, y1 = a*np.cos(T1), b*np.sin(T1)

plt.figure()
plt.scatter(x1, y1, c='k', s=5)
plt.gca().set_aspect('equal')
plt.show()

dx1 = x1[1:] - x1[:-1]
dy1 = y1[1:] - y1[:-1]
H1 = np.sqrt(dx1**2 + dy1**2)

plt.figure()
plt.axhline(y=h, c='gray', linewidth=1, zorder=1)
plt.plot(T1[1:], H1, c='k', zorder=2)
plt.show()

# 1/kappa sampled points on an ellipse

h = 0.1
N = int(p//h)
if N % 2 == 1: N += 1
T = np.linspace(0, 2*np.pi, N + 1)
x, y = a*np.cos(T), b*np.sin(T)
dx = x[1:] - x[:-1]
dy = y[1:] - y[:-1]
D = np.concatenate([[0], np.cumsum(np.sqrt(dx**2 + dy**2))])
nx, ny = -a*np.cos(T), -b*np.sin(T)
kappa = np.sqrt(nx**2 + ny**2)
rho = 1/kappa
S = np.concatenate([[0], np.cumsum(rho)])
S *= D.max()/S.max()
T2 = np.interp(S, D, T)
x2, y2 = a*np.cos(T2), b*np.sin(T2)

plt.figure()
plt.scatter(x2, y2, c='k', s=5)
plt.gca().set_aspect('equal')
plt.show()

dx2 = x2[1:] - x2[:-1]
dy2 = y2[1:] - y2[:-1]
H2 = np.sqrt(dx2**2 + dy2**2)

plt.figure()
plt.axhline(y=h, c='gray', linewidth=1, zorder=1)
plt.plot(T2[1:], H2, c='k', zorder=2)
plt.show()
