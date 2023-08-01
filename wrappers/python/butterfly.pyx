# cython: language_level=3, embedsignature=True, auto_pickle=False

import numpy as np

from enum import Enum

from bbox cimport *
from defs cimport *
from ellipse cimport *
from fac_helm2 cimport *
from geom cimport *
from helm2 cimport *
from layer_pot cimport *
from linalg cimport *
from mat cimport *
from mat_block_dense cimport *
from mat_dense_complex cimport *
from mat_diag_real cimport *
from mat_func cimport *
from mat_identity cimport *
from mat_product cimport *
from perm cimport *
from points cimport *
from ptr_array cimport *
from quadtree cimport *
from quadtree_node cimport *
from rand cimport *
from real_array cimport *
from size_array cimport *
from tree cimport *
from tree_node cimport *
from types cimport *
from vec cimport *
from vectors cimport *

def seed(BfSize seed):
    bfSeed(seed)

cdef reify_mat(BfMat *mat):
    cdef BfType type_ = bfMatGetType(mat)
    if type_ == BF_TYPE_MAT_BLOCK_DENSE:
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
    def xy(self):
        return (self.xmin, self.ymin)

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

    @property
    def dx(self):
        return self.xmax - self.xmin

    @property
    def dy(self):
        return self.ymax - self.ymin

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
    def make_multilevel(Helm2 helm, Quadtree srcTree, Quadtree tgtTree=None):
        return reify_mat(bfFacHelm2MakeMultilevel(
            &helm.helm,
            srcTree.quadtree,
            srcTree.quadtree if tgtTree is None else tgtTree.quadtree))

cdef class Helm2:
    cdef BfHelm2 helm

    def __init__(self, k, layerPot, alpha=None, beta=None):
        self.helm.k = k
        self.helm.layerPot = layerPot.value
        if layerPot == LayerPot.CombinedField:
            if alpha is None or beta is None:
                raise ValueError('alpha and beta must both be set if layerPot is CombinedField')
            self.helm.alpha.real = alpha.real
            self.helm.alpha.imag = alpha.imag
            self.helm.beta.real = beta.real
            self.helm.beta.imag = beta.imag
        else:
            if alpha is not None or beta is not None:
                raise ValueError('only pass alpha and beta when layerPot is CombinedField')

    @property
    def k(self):
        return self.helm.k

    @property
    def layer_pot(self):
        return self.helm.layerPot

    @property
    def alpha(self):
        return self.helm.alpha

    @property
    def beta(self):
        return self.helm.beta

    def get_kernel_matrix(self, Points2 Xsrc, Points2 Xtgt=None, Vectors2 Nsrc=None, Vectors2 Ntgt=None):
        return reify_mat(bfHelm2GetKernelMatrix(
            &self.helm,
            Xsrc.points,
            NULL if Xtgt is None else Xtgt.points,
            NULL if Nsrc is None else Nsrc.vectors,
            NULL if Ntgt is None else Ntgt.vectors))

    def apply_KR_correction(self, Mat mat, BfSize KR_order, Points2 X, Vectors2 N, Tree tree=None):
        if tree is None:
            bfHelm2ApplyKrCorrection(&self.helm, KR_order, X.points, N.vectors, mat.mat)
        else:
            bfHelm2ApplyKrCorrectionTree(&self.helm, KR_order, X.points, N.vectors, tree.tree, mat.mat)

    def apply_block_KR_correction(self, Mat mat, offsets, BfSize KR_order, Points2 X, Vectors2 N, Tree tree=None):
        cdef BfSizeArray *offsets_ = bfSizeArrayNewWithCapacity(len(offsets))
        for offset in offsets:
            bfSizeArrayAppend(offsets_, offset)
        if tree is None:
            bfHelm2ApplyBlockCorrection(&self.helm, offsets_, KR_order, X.points, N.vectors, mat.mat)
        else:
            bfHelm2ApplyBlockCorrectionTree(&self.helm, offsets_, KR_order, X.points, N.vectors, tree.tree, mat.mat)
        bfSizeArrayDeinitAndDealloc(&offsets_)

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

def solve_gmres(Mat A, Mat B, Mat X0=None, float tol=1e-15, max_num_iter=None, Mat M=None):
    cdef BfSize num_iter
    Y = reify_mat(bfSolveGMRES(
        A.mat,
        B.mat,
        NULL if X0 is None else X0.mat,
        tol,
        A.shape[0] - 1 if max_num_iter is None else max_num_iter,
        &num_iter,
        NULL if M is None else M.mat))
    return Y, num_iter

cdef class Mat:
    cdef BfMat *mat

    @property
    def shape(self):
        return (bfMatGetNumRows(self.mat), bfMatGetNumCols(self.mat))

    def __iadd__(self, Mat mat):
        bfMatAddInplace(self.mat, mat.mat)
        return self

    def to_mat_dense_complex(self):
        return MatDenseComplex.from_ptr(
            bfMatToMatDenseComplex(
                bfMatToType(self.mat, BF_TYPE_MAT_DENSE_COMPLEX)))

cdef class MatBlockDense(Mat):
    cdef BfMatBlockDense *mat_block_dense

    @staticmethod
    cdef MatBlockDense from_ptr(BfMatBlockDense *mat_block_dense):
        _ = MatBlockDense()
        _.mat = bfMatBlockDenseToMat(mat_block_dense)
        _.mat_block_dense = mat_block_dense
        return _

    def scale_cols(self, RealArray real_array):
        cdef BfVec *vec = bfRealArrayGetVecView(real_array.real_array)
        bfMatBlockDenseScaleCols(self.mat_block_dense, vec)
        bfVecDelete(&vec)

    def scale_cols(self, Vec vec):
        raise NotImplementedError()

    def get_blocks(self, I, J):
        cdef BfSize mBlk = len(I)
        cdef BfSize nBlk = len(J)
        cdef BfPtrArray blocks
        bfInitPtrArray(&blocks, mBlk*nBlk)
        cdef BfSize i
        cdef BfSize j
        cdef BfMat *block
        for i in range(mBlk):
            for j in range(nBlk):
                block = bfMatBlockDenseGetBlock(self.mat_block_dense, I[i], J[j])
                bfPtrArrayAppend(&blocks, block)
        blocks_ = MatBlockDense()
        blocks_.mat_block_dense = \
            bfMatBlockDenseNewFromBlocks(mBlk, nBlk, &blocks, BF_POLICY_VIEW)
        blocks_.mat = bfMatBlockDenseToMat(blocks_.mat_block_dense)
        return blocks_

cdef class MatDenseComplex(Mat):
    cdef BfMatDenseComplex *mat_dense_complex

    cdef size_t _buf_itemsize
    cdef Py_ssize_t[2] _buf_shape
    cdef Py_ssize_t[2] _buf_strides

    cdef void _buf_init(self):
        self._buf_itemsize = sizeof(BfComplex)
        self._buf_shape[0] = self.shape[0]
        self._buf_shape[1] = self.shape[1]
        self._buf_strides[1] = self._buf_itemsize
        self._buf_strides[0] = self.shape[1]*self._buf_strides[1]

    @staticmethod
    cdef MatDenseComplex from_ptr(BfMatDenseComplex *mat_dense_complex):
        _ = MatDenseComplex()
        _.mat = bfMatDenseComplexToMat(mat_dense_complex)
        _.mat_dense_complex = mat_dense_complex
        _._buf_init()
        return _

    def __getbuffer__(self, Py_buffer *buf, int flags):
        buf.buf = <char *>bfMatDenseComplexGetDataPtr(self.mat_dense_complex)
        buf.format = 'Zd'
        buf.internal = <void *>NULL
        buf.itemsize = self._buf_itemsize
        buf.len = self._buf_shape[0]*self._buf_shape[1]*self._buf_itemsize
        buf.ndim = 2
        buf.obj = self
        buf.readonly = 1
        buf.shape = self._buf_shape
        buf.strides = self._buf_strides
        buf.suboffsets = NULL

    def __releasebuffer__(self, Py_buffer *buf):
        pass

cdef class MatDiagReal(Mat):
    cdef BfMatDiagReal *matDiagReal

    @staticmethod
    cdef MatDiagReal from_constant(BfSize m, BfSize n, BfReal diag_value):
        mat_diag_real = MatDiagReal()
        mat_diag_real.matDiagReal = bfMatDiagRealNewConstant(m, n, diag_value)
        mat_diag_real.mat = bfMatDiagRealToMat(mat_diag_real.matDiagReal)
        return mat_diag_real

cdef class MatFunc(Mat):
    cdef BfMatFunc *MatFunc

cdef class MatIdentity(Mat):
    cdef BfMatIdentity *matIdentity

    def __cinit__(self, BfSize n):
        self.matIdentity = bfMatIdentityNew()
        self.mat = bfMatIdentityToMat(self.matIdentity)
        bfMatIdentityInit(self.matIdentity, n)

    def __truediv__(self, BfReal denom):
        m, n = self.shape
        return MatDiagReal.from_constant(m, n, 1/denom)

cdef class MatProduct(Mat):
    cdef BfMatProduct *matProduct

    def __cinit__(self):
        self.matProduct = bfMatProductNew()
        self.mat = bfMatProductToMat(self.matProduct)
        bfMatProductInit(self.matProduct)

    @staticmethod
    def from_factors(*factors):
        matProduct = MatProduct()
        for factor in factors:
            matProduct.post_multiply(factor)
        return matProduct

    cdef post_multiply(self, Mat mat):
        bfMatProductPostMultiply(self.matProduct, mat.mat)

cdef class Perm:
    cdef BfPerm *perm

    def __getitem__(self, i):
        return bfPermGetIndex(self.perm, i)

    def get_reverse(self):
        perm = Perm()
        perm.perm = bfPermGetReversePerm(self.perm)
        return perm

cdef class Points2:
    cdef BfPoints2 *points

    cdef Py_ssize_t[2] shape
    cdef Py_ssize_t[2] strides

    def __init__(self, BfSize capacity=16):
        self.points = bfPoints2NewWithCapacity(capacity)

    @staticmethod
    def from_point(point):
        if len(point) != 2:
            raise ValueError('len(point) != 2')
        cdef Points2 points = Points2.__new__(Points2)
        points.points = bfPoints2NewWithCapacity(1)
        cdef BfPoint2 point_ = [point[0], point[1]]
        bfPoints2Append(points.points, point_)
        return points

    @staticmethod
    def sample_poisson_disk(Bbox2 bbox, BfReal min_dist, BfSize k=30):
        cdef Points2 points = Points2.__new__(Points2)
        points.points = bfPoints2SamplePoissonDisk(&bbox.bbox, min_dist, k)
        return points

    def __getbuffer__(self, Py_buffer *buf, int flags):
        itemsize = sizeof(BfReal)
        num_points = len(self)

        self.shape[0] = num_points
        self.shape[1] = 2

        self.strides[1] = itemsize
        self.strides[0] = 2*self.strides[1]

        buf.buf = <char *>bfPoints2GetDataPtr(self.points)
        buf.format = 'd'
        buf.internal = NULL
        buf.itemsize = itemsize
        buf.len = self.shape[0]*self.shape[1]*itemsize
        buf.ndim = 2
        buf.obj = self
        buf.readonly = 1
        buf.shape = self.shape
        buf.strides = self.strides
        buf.suboffsets = NULL

    def __releasebuffer__(self, Py_buffer *buf):
        pass

    def __len__(self):
        return bfPoints2GetSize(self.points)

    def __getitem__(self, i):
        cdef BfPoint2 point
        if i > len(self):
            raise IndexError()
        bfPoints2Get(self.points, i, point)
        return Point2(point[0], point[1])

    def extend(self, Points2 points):
        bfPoints2Extend(self.points, points.points)

cdef class Quadtree(Tree):
    cdef BfQuadtree *quadtree

    def __init__(self):
        raise RuntimeError("use factory functions to instantiate Quadtree")

    @staticmethod
    cdef from_ptr(BfQuadtree *quadtree):
        cdef Quadtree _ = Quadtree.__new__(Quadtree)
        _.quadtree = quadtree
        _.tree = bfQuadtreeToTree(_.quadtree)
        return _

    @staticmethod
    def from_points_and_normals(Points2 points, Vectors2 normals):
        cdef Quadtree quadtree = Quadtree.__new__(Quadtree)
        quadtree.quadtree = bfQuadtreeNew()
        bfQuadtreeInit(quadtree.quadtree, points.points, normals.vectors)
        quadtree.tree = bfQuadtreeToTree(quadtree.quadtree)
        return quadtree

    def get_level_nodes(self, BfSize level):
        cdef BfPtrArray levelPtrArray = bfTreeGetLevelPtrArray(self.tree, level)
        nodes = []
        cdef BfQuadtreeNode *node = NULL
        for i in range(bfPtrArraySize(&levelPtrArray)):
            node = <BfQuadtreeNode *>bfPtrArrayGet(&levelPtrArray, i)
            node_ = QuadtreeNode.from_ptr(node)
            nodes.append(node_)
        bfPtrArrayDeinit(&levelPtrArray)
        return nodes

cdef class QuadtreeNode(TreeNode):
    cdef BfQuadtreeNode *quadtree_node

    @staticmethod
    cdef from_ptr(BfQuadtreeNode *quadtree_node):
        _ = QuadtreeNode()
        _.quadtree_node = quadtree_node
        _.tree_node = bfQuadtreeNodeToTreeNode(quadtree_node)
        return _

    @property
    def quadtree(self):
        cdef BfQuadtree *quadtree = bfQuadtreeNodeGetQuadtree(self.quadtree_node)
        return Quadtree.from_ptr(quadtree)

    @property
    def parent(self):
        cdef BfQuadtreeNode *parent = \
            bfTreeNodeToQuadtreeNode(bfTreeNodeGetParent(self.tree_node))
        return QuadtreeNode.from_ptr(parent)

    @property
    def split(self):
        cdef BfPoint2 split
        bfQuadtreeNodeGetSplit(self.quadtree_node, split)
        return (split[0], split[1])

    @property
    def bbox(self):
        cdef BfBbox2 bbox = bfQuadtreeNodeGetBbox(self.quadtree_node)
        return Bbox2(bbox.min[0], bbox.max[0], bbox.min[1], bbox.max[1])

    def get_inds(self):
        cdef BfSize i0 = self.get_first_index()
        cdef BfSize i1 = self.get_last_index()
        cdef SizeArray inds = SizeArray.__new__(SizeArray)
        inds.sizeArray = bfSizeArrayNewWithCapacity(i1 - i0)
        cdef Perm perm = self.quadtree.perm
        cdef BfSize i
        for i in range(i0, i1):
            bfSizeArrayAppend(inds.sizeArray, perm[i])
        return inds

    def get_points(self):
        cdef BfQuadtree *quadtree = bfQuadtreeNodeGetQuadtree(self.quadtree_node)
        points = Points2()
        points.points = bfQuadtreeNodeGetPoints(self.quadtree_node, quadtree)
        return points

def sample_uniform(lo=0, hi=1, n=1, dtype=np.float64):
    if dtype == np.float64:
        samples = [(hi - lo)*bfRealUniform1() + lo for _ in range(n)]
    else:
        raise NotImplementedError()
    return samples[0] if n == 1 else samples

cdef class RealArray:
    cdef BfRealArray *real_array

    cdef Py_ssize_t[1] _buf_shape
    cdef Py_ssize_t[1] _buf_strides

    def __cinit__(self):
        self.real_array = bfRealArrayNew()

    def __init__(self):
        bfRealArrayInitWithDefaultCapacity(self.real_array)

    def __len__(self):
        return bfRealArrayGetSize(self.real_array)

    def __getbuffer__(self, Py_buffer *buf, int flags):
        itemsize = sizeof(BfReal)
        n = len(self)
        self._buf_shape[0] = n
        self._buf_strides[0] = itemsize
        buf.buf = <char *>bfRealArrayGetDataPtr(self.real_array)
        buf.format = 'd'
        buf.internal = NULL
        buf.itemsize = itemsize
        buf.len = n*itemsize
        buf.ndim = 1
        buf.obj = self
        buf.readonly = 1
        buf.shape = self._buf_shape
        buf.strides = self._buf_strides
        buf.suboffsets = NULL

    def __releasebuffer__(self, Py_buffer *buf):
        pass

    def __getitem__(self, index):
        if isinstance(index, Perm):
            copy = self.copy()
            copy.permute(index)
            return copy
        else:
            raise NotImplementedError(f'type {type(index)} currently unsupported')

    def copy(self):
        cdef RealArray _ = RealArray.__new__(RealArray)
        _.real_array = bfRealArrayNew()
        bfRealArrayInitCopy(_.real_array, self.real_array)
        return _

    def extend(self, RealArray arr):
        bfRealArrayExtend(self.real_array, arr.real_array)

    def permute(self, Perm perm):
        bfRealArrayPermute(self.real_array, perm.perm)

cdef class SizeArray:
    cdef BfSizeArray *sizeArray

    cdef Py_ssize_t[1] _buf_shape
    cdef Py_ssize_t[1] _buf_strides

    @staticmethod
    cdef from_ptr(BfSizeArray *sizeArray):
        cdef SizeArray _ = SizeArray.__new__(SizeArray)
        _.sizeArray = sizeArray
        return _

    def __len__(self):
        return bfSizeArrayGetSize(self.sizeArray)

    def __getbuffer__(self, Py_buffer *buf, int flags):
        itemsize = sizeof(BfSize)
        n = len(self)
        self._buf_shape[0] = n
        self._buf_strides[0] = itemsize
        buf.buf = <char *>bfSizeArrayGetDataPtr(self.sizeArray)
        buf.format = 'Q'
        buf.internal = NULL
        buf.itemsize = itemsize
        buf.len = n*itemsize
        buf.ndim = 1
        buf.obj = self
        buf.readonly = 1
        buf.shape = self._buf_shape
        buf.strides = self._buf_strides
        buf.suboffsets = NULL

    def __releasebuffer__(self, Py_buffer *buf):
        pass

cdef class Tree:
    cdef BfTree *tree

    @property
    def perm(self):
        perm = Perm()
        perm.perm = bfPermGetView(bfTreeGetPerm(self.tree))
        return perm

cdef class TreeNode:
    cdef BfTreeNode *tree_node

    def get_num_points(self):
        return bfTreeNodeGetNumPoints(self.tree_node)

    def get_first_index(self):
        return bfTreeNodeGetFirstIndex(self.tree_node)

    def get_last_index(self):
        return bfTreeNodeGetLastIndex(self.tree_node)

cdef class Vec:
    cdef BfVec *vec

cdef class Vectors2:
    cdef BfVectors2 *vectors

    def __cinit__(self):
        self.vectors = bfVectors2NewEmpty()

    def __len__(self):
        return bfVectors2GetSize(self.vectors)

    def extend(self, Vectors2 vectors):
        bfVectors2Extend(self.vectors, vectors.vectors)
