cdef extern from "bf/types.h":
    cdef struct BfPtrArray:
        pass

    struct BfMat:
        pass

    struct BfMatBlockCoo:
        pass

    struct BfMatBlockDense:
        pass

    struct BfMatBlockDiag:
        pass

    struct BfMatCsrReal:
        pass

    struct BfMatDenseComplex:
        pass

    struct BfMatDenseReal:
        pass

    struct BfMatDiagReal:
        pass

    struct BfMatDiff:
        pass

    struct BfMatFunc:
        pass

    struct BfMatIdentity:
        pass

    struct BfMatProduct:
        pass

    struct BfMatPython:
        pass

    struct BfVec:
        pass

    struct BfTree:
        pass

    struct BfQuadtree:
        pass

    struct BfTreeNode:
        pass

    struct BfTreeLevelIter:
        pass

    struct BfNodeSpan:
        pass

    struct BfFacStreamer:
        pass

    struct BfFac:
        pass

    struct BfTrimesh:
        pass

    cdef enum BfTypes:
        # Mat hierarchy
        BF_TYPE_MAT
        BF_TYPE_MAT_COO_COMPLEX
        BF_TYPE_MAT_COO_REAL
        BF_TYPE_MAT_CSR_REAL
        BF_TYPE_MAT_DIAG_REAL
        BF_TYPE_MAT_DIFF
        BF_TYPE_MAT_FUNC
        BF_TYPE_MAT_GIVENS_COMPLEX
        BF_TYPE_MAT_IDENTITY
        BF_TYPE_MAT_PERM
        BF_TYPE_MAT_PYTHON
        BF_TYPE_MAT_PRODUCT
        BF_TYPE_MAT_SUM
        BF_TYPE_MAT_ZERO

        # MatBlock hierarchy:
        BF_TYPE_MAT_BLOCK
        BF_TYPE_MAT_BLOCK_COO
        BF_TYPE_MAT_BLOCK_DENSE
        BF_TYPE_MAT_BLOCK_DIAG

        # MatDense hierarchy:
        BF_TYPE_MAT_DENSE
        BF_TYPE_MAT_DENSE_COMPLEX
        BF_TYPE_MAT_DENSE_REAL

        # Lu hierarchy:
        BF_TYPE_LU
        BF_TYPE_LU_CSR_REAL
        BF_TYPE_LU_DENSE_COMPLEX

        # Vec hierarchy
        BF_TYPE_VEC
        BF_TYPE_VEC_COMPLEX
        BF_TYPE_VEC_REAL
        BF_TYPE_VEC_ZERO

        # Tree hierarchy
        BF_TYPE_TREE
        BF_TYPE_FIEDLER_TREE
        BF_TYPE_INTERVAL_TREE
        BF_TYPE_OCTREE
        BF_TYPE_QUADTREE

        # Tree node hierarchy
        BF_TYPE_TREE_NODE
        BF_TYPE_FIEDLER_TREE_NODE
        BF_TYPE_INTERVAL_TREE_NODE
        BF_TYPE_OCTREE_NODE
        BF_TYPE_QUADTREE_NODE

        # Tree iterator hierarchy
        BF_TYPE_TREE_ITER
        BF_TYPE_TREE_ITER_POST_ORDER

        # The total number of types
        BF_TYPE_COUNT

        # Placeholder for missing type information
        BF_TYPE_UNKNOWN

    ctypedef BfTypes BfType

type_to_str = {
    BF_TYPE_LU: 'Lu'
}
