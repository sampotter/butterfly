import colorcet as cc
import itertools as it
import matplotlib.pyplot as plt
import numpy as np
import scipy.special

from util import complex_to_hsv

plt.ion()

MAKE_DEBUG_PLOTS = True
MAKE_FIELD_PLOT = True

order = 10
n = 800
nplot = 50 # number of grid points along x and y axes in plot
nerr = 16 # number of points on circle at which to sample relative error

k = 30 # wavenumber
lam = 2*np.pi/k # wavelength

x0, y0 = 1, 0.2
q0 = np.array([x0, y0])
a, b = 1.5, 2
alpha = np.array([a, b])

T = np.linspace(0, 2*np.pi, n, endpoint=False)
h = 2*np.pi/n

x = a*np.cos(T) + x0
y = b*np.sin(T) + y0
Q = np.array([x, y]).T

dx, d2x = -a*np.sin(T), -a*np.cos(T)
dy, d2y =  b*np.cos(T), -b*np.sin(T)
dsdt = np.sqrt(dx**2 + dy**2)
tx, ty = dx/dsdt, dy/dsdt
nx, ny = ty, -tx # outward-facing
N = np.array([nx, ny]).T
L = np.cumsum(h*dsdt)

kappa = (dx*d2y - dy*d2x)/dsdt**3

def H0(x):
    return scipy.special.j0(x) + 1j*scipy.special.y0(x)

def dH0(x):
    return -scipy.special.j1(x) - 1j*scipy.special.y1(x)

def G(p, q, k):
    r = np.linalg.norm(p - q)
    return 1j*H0(k*r)/4

def gradG(p, q, k):
    r = np.linalg.norm(p - q)
    return k*1j*dH0(k*r)*(q - p)/(4*r)

nwaves = int(np.floor(L[-1]/lam))
print(f'number of wavelengths in curve: {nwaves}')

xsrc, ysrc = x0 + 0.1, y0 - 0.2
psrc = np.array([xsrc, ysrc])

# set up Neumann BCs
lhs = np.array([n@gradG(q, psrc, k) for n, q in zip(N, Q)])

# set up exact solution
def fexact(p):
    return G(p, psrc, k)

# Kapur-Rokhlin quadrature nodes
wKR = {2: [1.825748064736159e0,
           -1.325748064736159e0],
       6: [4.967362978287758e0,
           -16.20501504859126e0,
           25.85153761832639e0,
           -22.22599466791883e0,
           9.930104998037539e0,
           -1.817995878141594e0],
       10: [7.832432020568779e0,
            -4.565161670374749e1,
            1.452168846354677e2,
            -2.901348302886379e2,
            3.870862162579900e2,
            -3.523821383570681e2,
            2.172421547519342e2,
            -8.707796087382991e1,
            2.053584266072635e1,
            -2.166984103403823e0]}


# set up system matrix
A = np.zeros((n, n), dtype=np.complex128)
for i, j in it.product(range(n), repeat=2):
    if i == j:
        A[i, j] = 1/2
    else:
        A[i, j] = N[i]@gradG(Q[i], Q[j], k)*h*dsdt[j]

# apply KR correction
w = wKR[order]
for i in range(n):
    for p in range(order):
        di = p + 1
        A[i, (i + di) % n] *= 1 + w[p]
        A[i, (i - di) % n] *= 1 + w[p]

# compute righthand side
rhs = scipy.linalg.solve(A, lhs)

# function for evaluating solution using single-layer potential
def f(p):
    if np.sum(((p - q0)/alpha)**2) < 1:
        return np.nan
    else:
        return sum(G(p, q, k)*rhs[i]*dsdt[i]*h for i, q in enumerate(Q))

# compute relative error
rerr = 3*max(a, b)
theta = np.linspace(0, 2*np.pi, nerr, endpoint=False)
Perr = np.array([x0 + rerr*np.cos(theta), y0 + rerr*np.sin(theta)]).T
relerr = np.array([abs(f(p) - fexact(p))/abs(fexact(p)) for p in Perr])
relerr_mean = np.mean(relerr)
relerr_bar = abs(relerr - relerr_mean).max()
print(f'relative error mean: {relerr_mean} (+/- {relerr_bar})')

if MAKE_DEBUG_PLOTS:
    ext = 1.1*max(abs(x).max(), abs(y).max())

    plt.figure()
    plt.plot(x.tolist() + [x[0]], y.tolist() + [y[0]], c='k', zorder=1)
    plt.quiver(x, y, nx, ny)
    plt.scatter(xsrc, ysrc, c='r')
    plt.xlim(-ext, ext)
    plt.ylim(-ext, ext)
    plt.gca().set_aspect('equal')
    plt.show()

    plt.figure()
    plt.plot(T, abs(1/kappa), label=r'$1/|\kappa|$')
    plt.legend()
    plt.show()

    plt.figure()
    plt.plot(T, np.real(rhs))
    plt.plot(T, np.imag(rhs))
    plt.show()

    plt.figure()
    plt.imshow(np.log10(abs(np.real(A))))
    plt.gca().set_aspect('equal')
    plt.show()

if MAKE_FIELD_PLOT:
    ext = 1.1*abs(Perr).ravel().max()

    xgrid, ygrid = np.linspace(-ext, ext, nplot), np.linspace(-ext, ext, nplot)
    X, Y = np.meshgrid(xgrid, ygrid, indexing='xy')
    P = np.array([X.ravel(), Y.ravel()]).T
    F = np.array([f(p) for p in P]).reshape(X.shape)

    F_exact = np.array([fexact(p) for p in P]).reshape(X.shape)

    error = abs(F - F_exact)/abs(F_exact)

    vmax = abs(F_exact).max()
    vmin = -vmax
    extent = [-ext, ext, ext, -ext]
    cmap = cc.cm.CET_D1A

    plt.figure(figsize=(12, 3))
    plt.subplot(1, 3, 1)
    plt.imshow(np.real(F), extent=extent, vmin=vmin, vmax=vmax, cmap=cmap)
    plt.colorbar()
    plt.scatter(*Perr.T, s=10, c='k', marker='x', zorder=2)
    plt.subplot(1, 3, 2)
    plt.imshow(np.real(F_exact), extent=extent, vmin=vmin, vmax=vmax, cmap=cmap)
    plt.colorbar()
    plt.scatter(*Perr.T, s=10, c='k', marker='x', zorder=2)
    plt.subplot(1, 3, 3)
    plt.imshow(np.log10(np.maximum(1e-16, error)), extent=extent,
               cmap=cc.cm.rainbow, zorder=1)
    plt.colorbar()
    plt.scatter(*Perr.T, s=10, c='k', marker='x', zorder=2)
    plt.tight_layout()
    plt.show()

# checking out the Schur complement of A

m = 400
A11 = A[:m, :m]
A21 = A[m:, :m]
A12 = A[:m, m:]
A22 = A[m:, m:]
Arefl = A21@np.linalg.solve(A11, A12)
Asc = A22 - Arefl

plt.figure()
plt.imshow(complex_to_hsv(Arefl[25:375, 25:375]))
plt.show()

# I = [(0, 400)] * 5
I = [(0, 400), (0, 200), (0, 100), (0, 50),  (0, 25)]
J = [(0, 25),  (0, 50), (0, 100), (0, 200), (0, 400)]
# J = [(0, 25)] * 5
plt.figure()
for (i0, i1), (j0, j1) in zip(I, J):
    Aij = Arefl[i0:i1, j0:j1]
    S = np.linalg.svd(Aij, full_matrices=False)[1]
    plt.semilogy(S, label=rf'$A({i0}:{i1}, {j0}:{j1})$', linewidth=2, zorder=2)
S = np.linalg.svd(Arefl, full_matrices=False)[1]
plt.semilogy(S, linewidth=1, c='k', linestyle='--', zorder=1)
plt.axvline(x=nwaves)
plt.legend()
plt.show()
