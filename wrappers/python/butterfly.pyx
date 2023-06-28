# cython: language_level=3, embedsignature=True, auto_pickle=False

from bbox cimport *
from ellipse cimport *
from points cimport *
from quadtree cimport *
from rand cimport *
from real_array cimport *
from vectors cimport *

def seed(BfSize seed):
    bfSeed(seed)

cdef class Bbox2:
    cdef BfBbox2 bbox

    def __cinit__(self, xmin, xmax, ymin, ymax):
        self.bbox.min[0] = xmin
        self.bbox.max[0] = xmax
        self.bbox.min[1] = ymin
        self.bbox.max[1] = ymax

    @property
    def xmin(self):
        return self.bbox.min[0]

    @property
    def xmax(self):
        return self.bbox.max[0]

    @property
    def ymin(self):
        return self.bbox.min[1]

    @property
    def ymax(self):
        return self.bbox.max[1]

    def __str__(self):
        return f'[{self.xmin}, {self.xmax}] x [{self.ymin}, {self.ymax}]'

cdef class Ellipse:
    cdef BfEllipse ellipse

    def __init__(self, semi_major_axis, semi_minor_axis, center, theta):
        self.ellipse.semiMajorAxis = semi_major_axis
        self.ellipse.semiMinorAxis = semi_minor_axis
        self.ellipse.center[0] = center[0]
        self.ellipse.center[1] = center[1]
        self.ellipse.theta = theta

    @property
    def perimeter(self):
        return bfEllipseGetPerimeter(&self.ellipse)

    def sample_linspaced(self, num_points):
        X = Points2()
        T = Vectors2()
        N = Vectors2()
        W = RealArray()
        bfEllipseSampleLinspaced(
            &self.ellipse, num_points, X.points, T.vectors, N.vectors, W.real_array)
        return X, T, N, W

cdef class Point2:
    cdef BfPoint2 point

    def __init__(self, BfReal x, BfReal y):
        self.point[0] = x
        self.point[1] = y

    @property
    def x(self):
        return self.point[0]

    @x.setter
    def x(self, value):
        self.point[0] = value

    @property
    def y(self):
        return self.point[1]

    @y.setter
    def y(self, value):
        self.point[1] = value

    def __repr__(self):
        return f'({self.x}, {self.y})'

    def __getitem__(self, i):
        if i != 0 and i != 1:
            raise IndexError('index should be 0 or 1')
        return self.point[i]

    def __setitem__(self, i, value):
        if i != 0 and i != 1:
            raise IndexError('index should be 0 or 1')
        self.point[i] = value

cdef class Points2:
    cdef BfPoints2 *points

    def __init__(self):
        self.points = bfPoints2NewEmpty()

    @staticmethod
    def sample_poisson_disk(Bbox2 bbox, BfReal min_dist, BfSize k=30):
        points = Points2()
        bfPoints2DeinitAndDealloc(&points.points) # TODO: remove this...
        points.points = bfPoints2SamplePoissonDisk(&bbox.bbox, min_dist, k)
        return points

    def __len__(self):
        return bfPoints2GetSize(self.points)

    def __getitem__(self, i):
        cdef BfPoint2 point
        bfPoints2Get(self.points, i, point)
        return Point2(point[0], point[1])

    def extend(self, Points2 points):
        bfPoints2Extend(self.points, points.points)

cdef class Quadtree:
    cdef BfQuadtree *quadtree

    def __cinit__(self):
        self.quadtree = bfQuadtreeNew()

    def __init__(self, Points2 points, Vectors2 normals):
        bfQuadtreeInit(self.quadtree, points.points, normals.vectors)

    def __dealloc__(self):
        bfQuadtreeDeinitAndDealloc(&self.quadtree)

cdef class RealArray:
    cdef BfRealArray *real_array

    def __cinit__(self):
        self.real_array = bfRealArrayNew()

    def __init__(self):
        bfRealArrayInitWithDefaultCapacity(self.real_array)

    def extend(self, RealArray arr):
        bfRealArrayExtend(self.real_array, arr.real_array)

cdef class Vectors2:
    cdef BfVectors2 *vectors

    def __cinit__(self):
        self.vectors = bfVectors2NewEmpty()

    def extend(self, Vectors2 vectors):
        bfVectors2Extend(self.vectors, vectors.vectors)
