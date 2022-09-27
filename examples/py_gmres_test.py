import matplotlib.pyplot as plt

plt.ion()

import numpy as np
import scipy.linalg

# from scipy.linalg.blas import zrotg

def zrotg(a, b):
    if a == 0:
        c, s = 0, 1
    else:
        tmp = b/a
        c = 1/np.sqrt(1 + abs(tmp)**2)
        s = tmp*c
    return c, np.conj(s)

n = 512
A = np.fromfile("A.bin", dtype=np.complex128).reshape(n, n)
b = np.fromfile("phi.bin", dtype=np.complex128)
x0 = np.zeros_like(b)

# n = 128
# A = np.random.randn(n, n) + 1j*np.random.randn(n, n)
# b = np.random.randn(n)
# x0 = np.random.randn(n)

x_gt = np.linalg.solve(A, b)

tol = 1e-12
res = dict()

M = [n]

for m in M:
    print(f'm = {m}')

    r = b - A@x0
    beta = np.linalg.norm(r)

    res[m] = [beta]

    while res[m][-1] > tol*beta:
        s = np.zeros(m + 1, dtype=np.complex128)
        s[0] = np.linalg.norm(r)

        V = np.zeros((r.size, m + 1), dtype=A.dtype)
        V[:, 0] = r/np.linalg.norm(r)
        H = np.zeros((m + 1, m), dtype=np.complex128)
        g_c, g_s = [], []

        for i in range(m):
            w = A@V[:, i]
            for k in range(i + 1):
                H[k, i] = V[:, k].conj().T@w
                w -= H[k, i]*V[:, k]
            H[i + 1, i] = np.linalg.norm(w)
            V[:, i+1] = w/H[i + 1, i]

            for j, (cos, sin) in enumerate(zip(g_c, g_s)):
                R = np.eye(i + 2, dtype=np.complex128)
                R[j:j+2, j:j+2] = np.array([[cos, -sin], [sin, cos]])
                H[:i+2, i] = R.conj().T@H[:i+2, i]

            cos, sin = zrotg(*H[i:i+2, i]); cos = np.real(cos)
            g_c.append(cos)
            g_s.append(sin)

            R = np.eye(i + 2, dtype=np.complex128)
            R[i:i+2, i:i+2] = np.array([[cos, -sin], [sin, cos]])
            H[:i+2, i] = R.conj().T@H[:i+2, i]

            s[:i+2] = R.conj().T@s[:i+2]

            res[m].append(abs(s[i+1]))
            if res[m][-1] <= tol*beta:
                break

        import ipdb; ipdb.set_trace()

        y = scipy.linalg.solve_triangular(H[:i+1, :i+1], s[:i+1])
        x = x0 + V[:, :i+1]@y
        r = b - A@x

plt.figure()
plt.imshow(abs(H[:-1]), interpolation='none')
plt.show()

plt.figure()
plt.imshow(abs(V.T@V.conj()), interpolation='none')
plt.show()

plt.figure()
for m in M:
    plt.semilogy(res[m]/beta, label=f'{m}')
plt.legend()
plt.show()
