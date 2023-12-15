# cython: language_level=3, embedsignature=True, auto_pickle=False

import matplotlib.pyplot as plt
import numpy as np

from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

from pathlib import Path

cimport numpy as cnp

cnp.import_array()

from enum import Enum

from bf cimport bfInit

bfInit()

from bbox cimport *
from defs cimport *
from ellipse cimport *
from fac cimport *
from fac_helm2 cimport *
from fac_streamer cimport *
from geom cimport *
from helm2 cimport *
from indexed_mat cimport *
from layer_pot cimport *
from linalg cimport *
from mat cimport *
from mat_block_coo cimport *
from mat_block_dense cimport *
from mat_block_diag cimport *
from mat_csr_real cimport *
from mat_dense_complex cimport *
from mat_dense_real cimport *
from mat_diag_real cimport *
from mat_diff cimport *
from mat_func cimport *
from mat_identity cimport *
from mat_product cimport *
from mat_python cimport *
from node_span cimport *
from perm cimport *
from points cimport *
from ptr_array cimport *
from quadtree cimport *
from quadtree_node cimport *
from rand cimport *
from real_array cimport *
from size_array cimport *
from tree cimport *
from tree_level_iter cimport *
from tree_node cimport *
from tree_traversals cimport *
from trimesh cimport *
from types cimport *
from vec cimport *
from vectors cimport *

def seed(BfSize seed):
    bfSeed(seed)

cdef reify_mat(BfMat *mat):
    cdef BfType type_ = bfMatGetType(mat)
    if type_ == BF_TYPE_MAT:
        raise RuntimeError()
    elif type_ == BF_TYPE_MAT_DIFF:
        return MatDiff.from_ptr(bfMatToMatDiff(mat))
    elif type_ == BF_TYPE_MAT_PRODUCT:
        return MatProduct.from_ptr(bfMatToMatProduct(mat))
    elif type_ == BF_TYPE_MAT_PYTHON:
        return MatPython.from_ptr(bfMatToMatPython(mat))
    elif type_ == BF_TYPE_MAT_BLOCK_DENSE:
        return MatBlockDense.from_ptr(bfMatToMatBlockDense(mat))
    elif type_ == BF_TYPE_MAT_DENSE_COMPLEX:
        return MatDenseComplex.from_ptr(bfMatToMatDenseComplex(mat))
    else:
        raise TypeError(f'failed to reify BfMat: got {type_}')

cdef reify_tree(BfTree *tree):
    cdef BfType type_ = bfTreeGetType(tree)
    if type_ == BF_TYPE_QUADTREE:
        return Quadtree.from_ptr(bfTreeToQuadtree(tree))
    else:
        raise TypeError(f'failed to reify BfTree: got {type_}')

cdef reify_tree_node(BfTreeNode *treeNode):
    cdef BfType type_ = bfTreeNodeGetType(treeNode)
    if type_ == BF_TYPE_QUADTREE_NODE:
        return QuadtreeNode.from_ptr(bfTreeNodeToQuadtreeNode(treeNode))
    else:
        raise TypeError(f'failed to reify BfTreeNode: got {type_}')

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

class Policy(Enum):
    View = BF_POLICY_VIEW
    Copy = BF_POLICY_COPY
    Steal = BF_POLICY_STEAL

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

cdef class FacStreamer:
    cdef BfFacStreamer *facStreamer

    @staticmethod
    def from_trees(Tree rowTree, Tree colTree,
                   rowTreeInitDepth=None,
                   colTreeInitDepth=None, BfReal tol=1e-15,
                   BfSize minNumRows=20, BfSize minNumCols=20,
                   bint compareRelativeErrors=False):
        if rowTreeInitDepth is None:
            rowTreeInitDepth = 1
        if colTreeInitDepth is None:
            colTreeInitDepth = colTree.get_max_depth()

        cdef BfFacSpec facSpec
        facSpec.rowTree = rowTree.tree
        facSpec.colTree = colTree.tree
        facSpec.rowTreeInitDepth = rowTreeInitDepth
        facSpec.colTreeInitDepth = colTreeInitDepth
        facSpec.tol = tol
        facSpec.minNumRows = minNumRows
        facSpec.minNumCols = minNumCols
        facSpec.compareRelativeErrors = compareRelativeErrors

        cdef FacStreamer _ = FacStreamer.__new__(FacStreamer)
        _.facStreamer = bfFacStreamerNew()
        bfFacStreamerInit(_.facStreamer, &facSpec)
        return _

    def _feed_ndarray(self, cnp.ndarray arr):
        self._feed_Mat(Mat.from_ndarray(arr))

    def _feed_Mat(self, Mat mat):
        bfFacStreamerFeed(self.facStreamer, mat.mat)

    def feed(self, mat):
        if isinstance(mat, np.ndarray):
            self._feed_ndarray(mat)
        elif isinstance(mat, Mat):
            self._feed_mat(mat)
        else:
            raise NotImplementedError()

    def is_done(self):
        return bfFacStreamerIsDone(self.facStreamer)

    def toMat(self):
        cdef BfFac *fac = bfFacStreamerGetFac(self.facStreamer)
        cdef MatProduct _ = MatProduct.__new__(MatProduct)
        _.matProduct = bfFacGetMatProduct(fac, BF_POLICY_STEAL)
        bfFacDeinitAndDealloc(&fac)
        _.mat = bfMatProductToMat(_.matProduct)
        return _

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

    @staticmethod
    cdef from_ndarray(cnp.ndarray arr):
        if arr.dtype == np.float64:
            return MatDenseReal.from_ndarray(arr)
        elif arr.dtype == np.complex128:
            return MatDenseComplex.from_ndarray(arr)
        else:
            raise NotImplementedError()

    @property
    def shape(self):
        return (bfMatGetNumRows(self.mat), bfMatGetNumCols(self.mat))

    def __iadd__(Mat self, Mat mat):
        bfMatAddInplace(self.mat, mat.mat)
        return self

    def _matmul_ndarray(self, cnp.ndarray arr):
        return self._matmul_mat(Mat.from_ndarray(arr))

    def _matmul_mat(self, Mat mat):
        return reify_mat(bfMatMul(self.mat, mat.mat))

    def __matmul__(self, mat):
        if isinstance(mat, np.ndarray):
            return self._matmul_ndarray(mat)
        elif isinstance(mat, Mat):
            return self._matmul_mat(mat)
        else:
            raise NotImplementedError()

    def _rmatmul_ndarray(self, cnp.ndarray arr):
        return self._rmatmul_mat(Mat.from_ndarray(arr))

    def _rmatmul_mat(self, Mat mat):
        return reify_mat(bfMatRmul(self.mat, mat.mat))

    def __rmatmul__(self, mat):
        if isinstance(mat, np.ndarray):
            return self._rmatmul_ndarray(mat)
        elif isinstance(mat, Mat):
            return self._rmatmul_mat(mat)
        else:
            raise NotImplementedError()

    def _sub_ndarray(self, cnp.ndarray arr):
        return self._sub_mat(Mat.from_ndarray(arr))

    def _sub_mat(self, Mat mat):
        return reify_mat(bfMatSub(self.mat, mat.mat))

    def __sub__(self, other):
        if isinstance(other, np.ndarray):
            return self._sub_ndarray(other)
        elif isinstance(other, Mat):
            return self._sub_mat(other)
        else:
            raise NotImplementedError()

    def to_mat_dense_complex(self):
        return MatDenseComplex.from_ptr(
            bfMatToMatDenseComplex(
                bfMatToType(self.mat, BF_TYPE_MAT_DENSE_COMPLEX)))

    def transpose(self):
        bfMatTranspose(self.mat)

cdef class MatBlockCoo(Mat):
    cdef BfMatBlockCoo *matBlockCoo

    @staticmethod
    def from_indexed_blocks(shape, indexedBlocks, policy):
        cdef BfSize numRows, numCols
        numRows, numCols = shape

        cdef BfPtrArray indexedBlocksPtrArray
        bfInitPtrArray(&indexedBlocksPtrArray, len(indexedBlocks))

        cdef Mat block
        cdef BfSize i0, j0
        cdef BfIndexedMat *indexedMat
        for i0, j0, arg in indexedBlocks:
            if isinstance(arg, np.ndarray):
                block = Mat.from_ndarray(arg)
            else:
                raise NotImplementedError()
            indexedMat = bfIndexedMatAlloc()
            indexedMat.i0 = i0
            indexedMat.j0 = j0
            indexedMat.mat = block.mat
            bfPtrArrayAppend(&indexedBlocksPtrArray, indexedMat)

        cdef MatBlockCoo _ = MatBlockCoo.__new__(MatBlockCoo)
        _.matBlockCoo = bfMatBlockCooNewFromIndexedBlocks(numRows, numCols, &indexedBlocksPtrArray, policy.value)
        _.mat = bfMatBlockCooToMat(_.matBlockCoo)

        return _

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

    def get_blocks(self, I, J, policy=Policy.View):
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
            bfMatBlockDenseNewFromBlocks(mBlk, nBlk, &blocks, policy.value)
        blocks_.mat = bfMatBlockDenseToMat(blocks_.mat_block_dense)
        return blocks_

    def get_row_blocks(self, I, policy=Policy.View):
        cdef BfSize numColBlocks = bfMatBlockDenseGetNumColBlocks(self.mat_block_dense)

        cdef BfPtrArray blocks
        bfInitPtrArray(&blocks, len(I)*numColBlocks)

        cdef BfMat *block = NULL
        cdef BfSize j
        for i in I:
            for j in range(numColBlocks):
                block = bfMatBlockDenseGetBlock(self.mat_block_dense, i, j)
                bfPtrArrayAppend(&blocks, block)

        cdef MatBlockDense _ = MatBlockDense.__new__(MatBlockDense)
        _.mat_block_dense = bfMatBlockDenseNewFromBlocks(len(I), numColBlocks, &blocks, policy.value)
        _.mat = bfMatBlockDenseToMat(_.mat_block_dense)

        return _

    def get_col_blocks(self, J, policy=Policy.View):
        cdef BfSize numRowBlocks = bfMatBlockDenseGetNumRowBlocks(self.mat_block_dense)

        cdef BfPtrArray blocks
        bfInitPtrArray(&blocks, numRowBlocks*len(J))

        cdef BfMat *block = NULL
        cdef BfSize i
        for i in range(numRowBlocks):
            for j in J:
                block = bfMatBlockDenseGetBlock(self.mat_block_dense, i, j)
                bfPtrArrayAppend(&blocks, block)

        cdef MatBlockDense _ = MatBlockDense.__new__(MatBlockDense)
        _.mat_block_dense = bfMatBlockDenseNewFromBlocks(numRowBlocks, len(J), &blocks, policy.value)
        _.mat = bfMatBlockDenseToMat(_.mat_block_dense)

        return _

cdef class MatBlockDiag(Mat):
    cdef BfMatBlockDiag *matBlockDiag

    @staticmethod
    def from_blocks(blocks, policy):
        cdef MatBlockDiag _ = MatBlockDiag.__new__(MatBlockDiag)

        cdef BfPtrArray blocksPtrArray
        bfInitPtrArray(&blocksPtrArray, len(blocks))

        cdef Mat block
        for block in blocks:
            bfPtrArrayAppend(&blocksPtrArray, block.mat)

        _.matBlockDiag = bfMatBlockDiagNewFromBlocks(&blocksPtrArray, policy.value)
        _.mat = bfMatBlockDiagToMat(_.matBlockDiag)

        bfPtrArrayDeinit(&blocksPtrArray)

        return _

cdef class MatCsrReal(Mat):
    cdef BfMatCsrReal *mat_csr_real

    @staticmethod
    def new_view_factor_matrix_from_trimesh(Trimesh trimesh, BfSize[::1] I=None, BfSize[::1] J=None):
        if I is None:
            I = np.arange(trimesh.num_faces, dtype=np.uintp)
        cdef BfSizeArray *rowInds = bfSizeArrayNewView(len(I), &I[0])

        if J is None:
            J = np.arange(trimesh.num_faces, dtype=np.uintp)
        cdef BfSizeArray *colInds = bfSizeArrayNewView(len(J), &J[0])

        cdef MatCsrReal _ = MatCsrReal.__new__(MatCsrReal)
        _.mat_csr_real = bfMatCsrRealNewViewFactorMatrixFromTrimesh(trimesh.trimesh, rowInds, colInds)

        bfSizeArrayDeinitAndDealloc(&rowInds)
        bfSizeArrayDeinitAndDealloc(&colInds)

        return _

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

    @staticmethod
    cdef MatDenseComplex from_ndarray(cnp.ndarray arr):
        # Make sure arr is a packed, 2D, row-major array of complex doubles:
        assert arr.ndim == 2
        assert arr.flags.c_contiguous
        assert arr.itemsize == 16

        # cdef BfSize m = arr.shape[0]
        # cdef BfSize n = arr.shape[1]
        # cdef BfComplex[:, :] data = arr

        # cdef MatDenseComplex _ = MatDenseComplex.__new__(MatDenseComplex)
        # _.mat_dense_complex = bfMatDenseComplexNewViewFromPtr(m, n, &data[0, 0])
        # _.mat = bfMatDenseComplexToMat(_.mat_dense_complex)
        # return _

        cdef MatDenseComplex _ = MatDenseComplex.__new__(MatDenseComplex)
        _.mat_dense_complex = bfMatDenseComplexNewViewFromPyArray(arr)
        _.mat = bfMatDenseComplexToMat(_.mat_dense_complex)
        return _

    def __array__(self):
        cdef BfComplex *data = bfMatDenseComplexGetDataPtr(self.mat_dense_complex)
        return np.asarray(<BfComplex[:self.shape[0], :self.shape[1]]>data)

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

cdef class MatDenseReal(Mat):
    cdef BfMatDenseReal *matDenseReal

    @staticmethod
    cdef MatDenseReal from_ndarray(cnp.ndarray arr):
        # Make sure arr is a packed, 2D, row-major array of doubles:
        assert arr.ndim == 2
        assert arr.flags.c_contiguous
        assert arr.itemsize == 8

        cdef BfSize m = arr.shape[0]
        cdef BfSize n = arr.shape[1]
        cdef BfReal[:, :] data = arr

        cdef MatDenseReal _ = MatDenseReal.__new__(MatDenseReal)
        _.matDenseReal = bfMatDenseRealNewViewFromPtr(m, n, &data[0, 0], n, 1)
        _.mat = bfMatDenseRealToMat(_.matDenseReal)
        return _

cdef class MatDiagReal(Mat):
    cdef BfMatDiagReal *matDiagReal

    @staticmethod
    cdef MatDiagReal from_constant(BfSize m, BfSize n, BfReal diag_value):
        mat_diag_real = MatDiagReal()
        mat_diag_real.matDiagReal = bfMatDiagRealNewConstant(m, n, diag_value)
        mat_diag_real.mat = bfMatDiagRealToMat(mat_diag_real.matDiagReal)
        return mat_diag_real

cdef class MatDiff(Mat):
    cdef BfMatDiff *matDiff

    def __init__(self, Mat first, Mat second, policy=Policy.View):
        self.matDiff = bfMatDiffNew(first.mat, second.mat, policy.value)
        self.mat = bfMatDiffToMat(self.matDiff)

    @staticmethod
    cdef from_ptr(BfMatDiff *matDiff):
        cdef MatDiff _ = MatDiff.__new__(MatDiff)
        _.matDiff = matDiff
        _.mat = bfMatDiffToMat(_.matDiff)
        return _

    @property
    def first(self):
        return reify_mat(bfMatDiffGetFirst(self.matDiff))

    @property
    def second(self):
        return reify_mat(bfMatDiffGetSecond(self.matDiff))

    def get_blocks(self, I, J, policy=Policy.View):
        firstBlocks = self.first.get_blocks(I, J, policy)
        secondBlocks = self.second.get_blocks(I, J, policy)
        return MatDiff(firstBlocks, secondBlocks, policy)

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
    cdef list _factors

    def __cinit__(self):
        self.matProduct = bfMatProductNew()
        self.mat = bfMatProductToMat(self.matProduct)
        bfMatProductInit(self.matProduct)

    @staticmethod
    cdef from_ptr(BfMatProduct *matProduct):
        cdef MatProduct _ = MatProduct.__new__(MatProduct)
        _.matProduct = matProduct
        _.mat = bfMatProductToMat(_.matProduct)
        return _

    @staticmethod
    def from_factors(*factors):
        matProduct = MatProduct()
        for factor in factors:
            matProduct.post_multiply(factor)
        return matProduct

    @property
    def factors(self):
        if not self._factors:
            self._factors = [
                reify_mat(bfMatProductGetFactor(self.matProduct, i))
                for i in range(bfMatProductNumFactors(self.matProduct))]
        return self._factors

    def get_blocks(self, I, J, policy=Policy.View):
        if len(self.factors) == 0:
            raise RuntimeError()
        elif len(self.factors) == 1:
            newFactors = [self.factors[0].get_blocks(I, J, policy)]
        else:
            newFactors = [self.factors[0].get_row_blocks(I, policy)] \
                + self.factors[1:-1] + [self.factors[-1].get_col_blocks(J, policy)]
        return MatProduct.from_factors(*newFactors)

    cdef post_multiply(self, Mat mat):
        assert mat.mat != NULL
        bfMatProductPostMultiply(self.matProduct, mat.mat)

cdef class NodeSpan:
    cdef BfNodeSpan *nodeSpan

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
    def from_tree(Tree tree):
        cdef Quadtree _ = Quadtree.__new__(Quadtree)
        _.tree = tree.tree
        _.quadtree = bfTreeToQuadtree(_.tree)
        return _

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

    def plot_node_boxes(self, ax=None):
        if ax is None:
            ax = plt.gca()
        rects = []
        for node in self.nodes:
            bbox = node.bbox
            dx, dy = node.bbox.dx, node.bbox.dy
            rect = Rectangle(bbox.xy, dy, dx)
            rects.append(rect)
        pc = PatchCollection(rects, facecolor='none', edgecolor='k')
        ax.add_collection(pc)

cdef class QuadtreeNode(TreeNode):
    cdef BfQuadtreeNode *quadtreeNode

    @staticmethod
    cdef from_ptr(BfQuadtreeNode *quadtreeNode):
        _ = QuadtreeNode()
        _.quadtreeNode = quadtreeNode
        _.treeNode = bfQuadtreeNodeToTreeNode(quadtreeNode)
        return _

    @property
    def split(self):
        cdef BfPoint2 split
        bfQuadtreeNodeGetSplit(self.quadtreeNode, split)
        return (split[0], split[1])

    @property
    def bbox(self):
        cdef BfBbox2 bbox = bfQuadtreeNodeGetBbox(self.quadtreeNode)
        return Bbox2(bbox.min[0], bbox.max[0], bbox.min[1], bbox.max[1])

    def get_inds(self):
        cdef BfSize i0 = self.get_first_index()
        cdef BfSize i1 = self.get_last_index()
        cdef SizeArray inds = SizeArray.__new__(SizeArray)
        inds.sizeArray = bfSizeArrayNewWithCapacity(i1 - i0)
        cdef Perm perm = self.tree.perm
        cdef BfSize i
        for i in range(i0, i1):
            bfSizeArrayAppend(inds.sizeArray, perm[i])
        return inds

    def get_points(self):
        cdef BfQuadtree *quadtree = bfQuadtreeNodeGetQuadtree(self.quadtreeNode)
        points = Points2()
        points.points = bfQuadtreeNodeGetPoints(self.quadtreeNode, quadtree)
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

    @staticmethod
    cdef _from_node_span_NodeSpan(NodeSpan nodeSpan):
        cdef Tree tree = Tree.__new__(Tree)
        cdef Perm perm = Perm.__new__(Perm)
        tree.tree = bfTreeNewFromNodeSpan(nodeSpan.nodeSpan, &perm.perm)
        return tree, perm

    @staticmethod
    cdef _from_node_span_list(list nodes):
        if not all(isinstance(_, TreeNode) for _ in nodes):
            raise ValueError()
        cdef BfPtrArray ptrArray
        bfInitPtrArray(&ptrArray, len(nodes))
        cdef TreeNode node
        for node in nodes:
            bfPtrArrayAppend(&ptrArray, node.treeNode)
        cdef BfNodeSpan *nodeSpan = bfNodeSpanNewFromPtrArray(&ptrArray, BF_POLICY_COPY)
        cdef Tree tree = Tree.__new__(Tree)
        cdef Perm perm = Perm.__new__(Perm)
        tree.tree = bfTreeNewFromNodeSpan(nodeSpan, &perm.perm)
        return tree, perm

    @staticmethod
    def from_node_span(nodeSpan):
        if isinstance(nodeSpan, NodeSpan):
            return Tree._from_node_span_NodeSpan(nodeSpan)
        elif isinstance(nodeSpan, list):
            return Tree._from_node_span_list(nodeSpan)
        else:
            raise NotImplementedError()

    @staticmethod
    def from_node(TreeNode node):
        cdef Tree tree = Tree.__new__(Tree)
        cdef Perm perm = Perm.__new__(Perm)
        tree.tree = bfTreeNewFromNode(node.treeNode, &perm.perm)
        return tree, perm

    @staticmethod
    def for_middle_fac(Tree templateTree, BfSize p):
        cdef Tree tree = Tree.__new__(Tree)
        tree.tree = bfTreeNewForMiddleFac(templateTree.tree, p)
        return tree

    @property
    def perm(self):
        perm = Perm()
        perm.perm = bfPermGetView(bfTreeGetPerm(self.tree))
        return perm

    @property
    def root(self):
        return reify_tree_node(<BfTreeNode *>bfTreeGetRootNode(self.tree))

    @property
    def nodes(self):
        return TreeLevelIter.from_tree(self)

    def get_level_nodes(self, BfSize level):
        cdef BfPtrArray *levelPtrArray = bfTreeGetLevelPtrArray(self.tree, level)
        nodes = [reify_tree_node(<BfTreeNode *>bfPtrArrayGet(levelPtrArray, i))
                 for i in range(bfPtrArraySize(levelPtrArray))]
        bfPtrArrayDeinitAndDealloc(&levelPtrArray)
        return nodes

    def __eq__(self, Tree other):
        return self.tree == other.tree

    def get_max_depth(self):
        return bfTreeGetMaxDepth(self.tree)

cdef class TreeLevelIter:
    cdef BfTreeLevelIter *treeLevelIter

    @staticmethod
    def from_tree(Tree tree):
        cdef TreeLevelIter _ = TreeLevelIter.__new__(TreeLevelIter)
        _.treeLevelIter = bfTreeLevelIterNewFromTree(BF_TREE_TRAVERSAL_UNKNOWN, tree.tree)
        return _

    def __iter__(self):
        cdef BfPtrArray *levelNodes = NULL
        cdef BfType treeNodeType
        while not bfTreeLevelIterIsDone(self.treeLevelIter):
            levelNodes = bfTreeLevelIterGetLevelNodes(self.treeLevelIter)
            for i in range(bfPtrArraySize(levelNodes)):
                yield reify_tree_node(<BfTreeNode *>bfPtrArrayGet(levelNodes, i))
            bfTreeLevelIterNext(self.treeLevelIter)

cdef class TreeNode:
    cdef BfTreeNode *treeNode

    def __hash__(self):
        return <Py_hash_t>self.treeNode

    def __eq__(self, TreeNode other):
        return self.treeNode == other.treeNode

    def get_num_points(self):
        return bfTreeNodeGetNumPoints(self.treeNode)

    def get_first_index(self):
        return bfTreeNodeGetFirstIndex(self.treeNode)

    def get_last_index(self):
        return bfTreeNodeGetLastIndex(self.treeNode)

    @property
    def tree(self):
        return reify_tree(bfTreeNodeGetTree(self.treeNode))

    @property
    def parent(self):
        return reify_tree_node(bfTreeNodeGetParent(self.treeNode))

    @property
    def children(self):
        return [
            reify_tree_node(bfTreeNodeGetChild(self.treeNode, i))
            for i in range(self.max_num_children)
            if bfTreeNodeHasChild(self.treeNode, i)
        ]

    @property
    def depth(self):
        return bfTreeNodeGetDepth(self.treeNode)

    @property
    def max_num_children(self):
        return bfTreeNodeGetMaxNumChildren(self.treeNode)

    @property
    def i0(self):
        return self.get_first_index()

    @property
    def i1(self):
        return self.get_last_index()

cdef class Trimesh:
    cdef BfTrimesh *trimesh

    @staticmethod
    def from_obj(path):
        path_byte_string = str(path).encode('UTF-8')
        cdef char *path_c_string = path_byte_string
        cdef Trimesh _ = Trimesh.__new__(Trimesh)
        _.trimesh = bfTrimeshNewFromObjFile(path_c_string)
        return _

    @property
    def num_faces(self):
        return bfTrimeshGetNumFaces(self.trimesh)

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

############################################################################
# Other stuff...

# - need to set mat for DenseLu
# - need to call bfMatInit and pass BfMatVtable
# - need to set MatMul entry of BfMatMul with a C or Cython function
#   that will look up and call DenseLu._Mul

cdef class MatPython(Mat):
    cdef BfMatPython *matPython

    def __init__(self, m, n):
        self.matPython = bfMatPythonNewFromPyObject(self, m, n)
        self.mat = bfMatPythonToMat(self.matPython)

    @staticmethod
    cdef from_ptr(BfMatPython *matPython):
        cdef MatPython _ = MatPython.__new__(MatPython)
        _.matPython = matPython
        _.mat = bfMatPythonToMat(_.matPython)
        return _
