import matplotlib.pyplot as plt
import matplotlib.patches
import numpy as np
import scipy.optimize

plt.ion()

xmin, xmax = 0, 100
ymin, ymax = 0, 100

rmin, rmax = 1, 3
dmin = 2

class Ellipse:
    def __init__(self, x0, y0, rx, ry, theta):
        self.x0 = x0
        self.y0 = y0
        self.rx = rx
        self.ry = ry
        self.theta = theta
        self.R = np.array([[np.cos(theta), np.sin(theta)],
                           [-np.sin(theta), np.cos(theta)]])

    def f(self, t):
        return self.R.T@np.array([
            self.rx*np.cos(t), self.ry*np.sin(t)]) + np.array([self.x0, self.y0])

    def dfdt(self, t):
        return self.R@np.array([-self.rx*np.sin(t), self.ry*np.cos(t)])

    def plot(self, ax=None, **kwargs):
        if ax is None:
            ax = plt.gca()
        ax.add_patch(matplotlib.patches.Ellipse(
            (self.x0, self.y0), 2*self.rx, 2*self.ry, np.rad2deg(self.theta),
            **kwargs))

    def get_bounding_box(self):
        xmin = scipy.optimize.minimize_scalar(lambda t: self.f(t)[0], bracket=[-np.pi/2, 5*np.pi/2]).fun
        xmax = -scipy.optimize.minimize_scalar(lambda t: -self.f(t)[0], bracket=[-np.pi/2, 5*np.pi/2]).fun
        ymin = scipy.optimize.minimize_scalar(lambda t: self.f(t)[1], bracket=[-np.pi/2, 5*np.pi/2]).fun
        ymax = -scipy.optimize.minimize_scalar(lambda t: -self.f(t)[1], bracket=[-np.pi/2, 5*np.pi/2]).fun
        return xmin, xmax, ymin, ymax

    def plot_bounding_box(self, ax=None, **kwargs):
        if ax is None:
            ax = plt.gca()
        xmin, xmax, ymin, ymax = self.get_bounding_box()
        xm, ym = xmin, ymin
        w, h = xmax - xmin, ymax - ymin
        ax.add_patch(matplotlib.patches.Rectangle((xm, ym), w, h, **kwargs))

    def point_is_in_bounding_box(self, x, y):
        xmin, xmax, ymin, ymax = self.get_bounding_box()
        return xmin <= x <= xmax and ymin <= y <= ymax

    def overlaps(self, other):
        def check_overlap(e1, e2):
            xmin, xmax, ymin, ymax = e1.get_bounding_box()
            return e2.point_is_in_bounding_box(xmin, ymin) \
                or e2.point_is_in_bounding_box(xmin, ymax) \
                or e2.point_is_in_bounding_box(xmax, ymin) \
                or e2.point_is_in_bounding_box(xmax, ymax)
        return check_overlap(self, other) and check_overlap(other, self)

    def get_min_dist(self, other):
        def p(st):
            s, t = st
            return self.f(s) - other.f(t)
        def fun(st):
            return p(st)@p(st)
        def jac(st):
            s, t = st
            return np.array([
                self.dfdt(s)@p(st)/fun(st),
                -other.dfdt(t)@p(st)/fun(st)])
        return np.sqrt(scipy.optimize.minimize(
            fun, (0, 0), jac=jac,
            bounds=[(-np.pi/2, 5*np.pi/2), (-np.pi/2, 5*np.pi/2)]).fun)

def get_random_ellipse():
    x0 = np.random.uniform(xmin, xmax)
    y0 = np.random.uniform(ymin, ymax)
    rx, ry = np.random.uniform(rmin, rmax, 2)
    theta = np.random.uniform(0, 2*np.pi)
    return Ellipse(x0, y0, rx, ry, theta)

# ellipses = [get_random_ellipse() for _ in range(100)]

ellipses = [get_random_ellipse()]
while len(ellipses) < 100:
    ellipse = get_random_ellipse()
    if any(ellipse.overlaps(_) for _ in ellipses):
        continue
    min_dist = min(ellipse.get_min_dist(_) for _ in ellipses)
    if min_dist <= dmin:
        continue
    ellipses.append(ellipse)

plt.figure(figsize=(8, 8))
for _ in ellipses:
    _.plot(edgecolor='k', facecolor='orange')
    _.plot_bounding_box(edgecolor='k', facecolor='none', linestyle='--')
plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)
plt.gca().set_aspect('equal')
plt.show()
