import colorcet as cc
import matplotlib.pyplot as plt
import numpy as np

from numpy.polynomial import Chebyshev
from scipy.optimize import brentq
from scipy.special import j0, y0

k = 10
eps = 1e-5

rmax = 10
R = np.linspace(0, rmax, 1000)

p = 4 # number of polys for each interval

nodes = []
m = 1
while True:
    node = np.pi*(m - 1/4)/k
    if node > rmax:
        break
    bracket = (node - np.pi/(2*k), node + np.pi/(2*k))
    node = brentq(lambda r: j0(k*r), *bracket, xtol=1e-15, rtol=1e-15)
    nodes.append(node)
    m += 1
nodes.insert(0, 0) # need to include zero in here as an additional bracket
if nodes[-1] != rmax:
    nodes.append(rmax)
nodes = np.array(nodes)

cheb = []
for m, (r0, r1) in enumerate(zip(nodes[:-1], nodes[1:])):
    get_cheb = lambda deg: Chebyshev.interpolate(lambda r: j0(k*r), deg, domain=[r0, r1])
    deg = 1
    cheb.append(get_cheb(deg))
    chebNext = get_cheb(deg + 1)
    while abs((cheb[-1] - chebNext).coef[-1]) > eps:
        cheb[-1] = chebNext
        deg += 1
        chebNext = get_cheb(deg + 1)

@np.vectorize
def get_m(r):
    m = int(np.floor(k*r/np.pi + 1/4))
    if r < nodes[m]:
        m -= 1
    assert nodes[m] <= r <= nodes[m + 1]
    return m

@np.vectorize
def j0_cheb(r):
    return cheb[get_m(r)](r)

M = np.array([get_m(r) for r in R])

plt.figure()
plt.axhline(y=0, c='#cccccc', linestyle='--', zorder=0)
plt.plot(R, j0(k*R), c='k', zorder=1)
plt.scatter(nodes, j0(k*nodes), s=10, c='k', zorder=1)
plt.scatter(R, M/300, c=cc.cm.glasbey(M/300), zorder=2, s=2)
plt.xlim(0, rmax)
plt.show()

plt.figure()
plt.semilogy(R, abs(j0(k*R) - j0_cheb(R)))
plt.show()
