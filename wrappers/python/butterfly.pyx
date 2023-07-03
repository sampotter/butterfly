# cython: language_level=3, embedsignature=True, auto_pickle=False

from enum import Enum

from bbox cimport *
from ellipse cimport *
from fac_helm2 cimport *
from geom cimport *
from helm2 cimport *
from layer_pot cimport *
from mat cimport *
from mat_block_dense cimport *
from mat_dense_complex cimport *
from perm cimport *
from points cimport *
from quadtree cimport *
from rand cimport *
from real_array cimport *
from tree cimport *
from types cimport *
from vectors cimport *

def seed(BfSize seed):
    bfSeed(seed)

cdef reify_mat(BfMat *mat):
    cdef BfType type_ = bfMatGetType(mat)
    if type_ == BF_TYPE_MAT:
        return Mat.from_ptr(mat)
    elif type_ == BF_TYPE_MAT_BLOCK_DENSE:
        return MatBlockDense.from_ptr(bfMatToMatBlockDense(mat))
    elif type_ == BF_TYPE_MAT_DENSE_COMPLEX:
        return MatDenseComplex.from_ptr(bfMatToMatDenseComplex(mat))
    else:
        raise TypeError(f'failed to reify BfMat: got {bfMatGetType(mat)}')

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

cdef class FacHelm2:
    def __init__(self):
        raise TypeError("FacHelm2 can't be instantiated")

    @staticmethod
    def make_multilevel(Quadtree srcTree, BfReal k, layerPot, Quadtree
                        tgtTree=None, alpha=None, beta=None):
        cdef const BfQuadtree *srcTree_ = srcTree.quadtree

        cdef const BfQuadtree *tgtTree_ = srcTree_ if tgtTree is None else tgtTree.quadtree

        cdef BfComplex alpha_
        if alpha is not None:
            alpha_.real = alpha.real
            alpha_.imag = alpha.imag

        cdef BfComplex beta_
        if beta is not None:
            beta_.real = beta.real
            beta_.imag = beta.imag

        cdef BfMat *mat = bfFacHelm2MakeMultilevel(
            srcTree_, tgtTree_, k, layerPot.value,
            NULL if alpha is None else &alpha_,
            NULL if beta is None else &beta_)

        return reify_mat(mat)

cdef class Helm2:
    def __init__(self):
        raise TypeError("Helm2 can't be instantiated")

    @staticmethod
    def get_kernel_matrix(Points2 Xsrc, BfReal k, layerPot, Points2
                          Xtgt=None, Vectors2 Nsrc=None, Vectors2
                          Ntgt=None, alpha=None, beta=None):
        cdef BfComplex alpha_
        if alpha is not None:
            alpha_.real = alpha.real
            alpha_.imag = alpha.imag

        cdef BfComplex beta_
        if beta is not None:
            beta_.real = beta.real
            beta_.imag = beta.imag

        cdef BfMat *mat = bfHelm2GetKernelMatrix(
            Xsrc.points,
            NULL if Xtgt is None else Xtgt.points,
            NULL if Nsrc is None else Nsrc.vectors,
            NULL if Ntgt is None else Ntgt.vectors,
            k,
            layerPot.value,
            NULL if alpha is None else &alpha_,
            NULL if beta is None else &beta_)

        return reify_mat(mat)

class LayerPot(Enum):
    Unknown = BF_LAYER_POTENTIAL_UNKNOWN
    Single = BF_LAYER_POTENTIAL_SINGLE
    PvDouble = BF_LAYER_POTENTIAL_PV_DOUBLE
    PvNormalDerivSingle = BF_LAYER_POTENTIAL_PV_NORMAL_DERIV_SINGLE
    PvNormalDerivDouble = BF_LAYER_POTENTIAL_PV_NORMAL_DERIV_DOUBLE
    CombinedField = BF_LAYER_POTENTIAL_COMBINED_FIELD

    # Aliases:
    S = Single
    Sp = PvNormalDerivSingle
    D = PvDouble
    Dp = PvNormalDerivDouble

cdef class Mat:
    cdef BfMat *mat

    @staticmethod
    cdef Mat from_ptr(BfMat *mat):
        _ = Mat()
        _.mat = mat
        return _

    def apply_KR_correction(self):
        print('hi!')

cdef class MatBlockDense(Mat):
    cdef BfMatBlockDense *mat_block_dense

    @staticmethod
    cdef MatBlockDense from_ptr(BfMatBlockDense *mat_block_dense):
        _ = MatBlockDense()
        _.mat = bfMatBlockDenseToMat(mat_block_dense)
        _.mat_block_dense = mat_block_dense
        return _

cdef class MatDenseComplex(Mat):
    cdef BfMatDenseComplex *mat_dense_complex

    @staticmethod
    cdef MatDenseComplex from_ptr(BfMatDenseComplex *mat_dense_complex):
        _ = MatDenseComplex()
        _.mat = bfMatDenseComplexToMat(mat_dense_complex)
        _.mat_dense_complex = mat_dense_complex
        return _

cdef class Perm:
    cdef BfPerm *perm

    def get_reverse(self):
        perm = Perm()
        perm.perm = bfPermGetReversePerm(self.perm)
        return perm

cdef class Points2:
    cdef BfPoints2 *points

    def __init__(self):
        self.points = bfPoints2NewEmpty()

    @staticmethod
    def from_point(point):
        if len(point) != 2:
            raise ValueError('len(point) != 2')
        points = Points2()
        cdef BfPoint2 point_ = [point[0], point[1]]
        bfPoints2Append(points.points, point_)
        return points

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

cdef class Quadtree(Tree):
    cdef BfQuadtree *quadtree

    def __cinit__(self):
        self.quadtree = bfQuadtreeNew()
        self.tree = bfQuadtreeToTree(self.quadtree)

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

cdef class Tree:
    cdef BfTree *tree

    @property
    def perm(self):
        perm = Perm()
        perm.perm = bfPermGetView(bfTreeGetPerm(self.tree))
        return perm

cdef class Vectors2:
    cdef BfVectors2 *vectors

    def __cinit__(self):
        self.vectors = bfVectors2NewEmpty()

    def extend(self, Vectors2 vectors):
        bfVectors2Extend(self.vectors, vectors.vectors)
