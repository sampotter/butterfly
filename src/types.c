#include <bf/types.h>

bool bfTypeDerivedFrom(BfType derived, BfType parent) {
  /* For each entry in this lookup table, we should have:
   *
   *   `_[SUBTYPE][TYPE] == true` if `SUBTYPE` is a subtype of `TYPE`
   *
   * When a new type is added, this lookup table needs to be updated.
   *
   * NOTE: static arrays are zero-initialized by default in C99. */
  static bool _[BF_TYPE_COUNT][BF_TYPE_COUNT] = {
    [BF_TYPE_MAT] = {
      [BF_TYPE_MAT] = true,
    },
    [BF_TYPE_MAT_COO_COMPLEX] = {
      [BF_TYPE_MAT] =             true,
      [BF_TYPE_MAT_COO_COMPLEX] = true
    },
    [BF_TYPE_MAT_COO_REAL] = {
      [BF_TYPE_MAT] =          true,
      [BF_TYPE_MAT_COO_REAL] = true
    },
    [BF_TYPE_MAT_CSR_REAL] = {
      [BF_TYPE_MAT] =          true,
      [BF_TYPE_MAT_CSR_REAL] = true
    },
    [BF_TYPE_MAT_DIAG_REAL] = {
      [BF_TYPE_MAT]           = true,
      [BF_TYPE_MAT_DIAG_REAL] = true,
    },
    [BF_TYPE_MAT_GIVENS_COMPLEX] = {
      [BF_TYPE_MAT]                = true,
      [BF_TYPE_MAT_GIVENS_COMPLEX] = true,
    },
    [BF_TYPE_MAT_PRODUCT] = {
      [BF_TYPE_MAT]         = true,
      [BF_TYPE_MAT_PRODUCT] = true
    },
    [BF_TYPE_MAT_SUM] = {
      [BF_TYPE_MAT]     = true,
      [BF_TYPE_MAT_SUM] = true
    },
    [BF_TYPE_MAT_ZERO] = {
      [BF_TYPE_MAT]      = true,
      [BF_TYPE_MAT_ZERO] = true
    },

    /** MatBlock hierarchy: */
    [BF_TYPE_MAT_BLOCK] = {
      [BF_TYPE_MAT]       = true,
      [BF_TYPE_MAT_BLOCK] = true,
    },
    [BF_TYPE_MAT_BLOCK_COO] = {
      [BF_TYPE_MAT]           = true,
      [BF_TYPE_MAT_BLOCK]     = true,
      [BF_TYPE_MAT_BLOCK_COO] = true,
    },
    [BF_TYPE_MAT_BLOCK_DENSE] = {
      [BF_TYPE_MAT]             = true,
      [BF_TYPE_MAT_BLOCK]       = true,
      [BF_TYPE_MAT_BLOCK_DENSE] = true,
    },
    [BF_TYPE_MAT_BLOCK_DIAG] = {
      [BF_TYPE_MAT]            = true,
      [BF_TYPE_MAT_BLOCK]      = true,
      [BF_TYPE_MAT_BLOCK_DIAG] = true,
    },

    /** MatDense hierarchy: */
    [BF_TYPE_MAT_DENSE] = {
      [BF_TYPE_MAT]       = true,
      [BF_TYPE_MAT_DENSE] = true
    },
    [BF_TYPE_MAT_DENSE_COMPLEX] = {
      [BF_TYPE_MAT]               = true,
      [BF_TYPE_MAT_DENSE]         = true,
      [BF_TYPE_MAT_DENSE_COMPLEX] = true,
    },
    [BF_TYPE_MAT_DENSE_REAL] = {
      [BF_TYPE_MAT]            = true,
      [BF_TYPE_MAT_DENSE]      = true,
      [BF_TYPE_MAT_DENSE_REAL] = true,
    },

    /** Vec hierarchy: */
    [BF_TYPE_VEC] = {
      [BF_TYPE_VEC] = true
    },
    [BF_TYPE_VEC_COMPLEX] = {
      [BF_TYPE_VEC] =         true,
      [BF_TYPE_VEC_COMPLEX] = true
    },
    [BF_TYPE_VEC_REAL] = {
      [BF_TYPE_VEC] =      true,
      [BF_TYPE_VEC_REAL] = true
    },
    [BF_TYPE_VEC_ZERO] = {
      [BF_TYPE_VEC] =      true,
      [BF_TYPE_VEC_ZERO] = true
    },

    /** Tree hierarchy: */
    [BF_TYPE_TREE] = {
      [BF_TYPE_TREE] = true
    },
    [BF_TYPE_FIEDLER_TREE] = {
      [BF_TYPE_TREE] =         true,
      [BF_TYPE_FIEDLER_TREE] = true
    },
    [BF_TYPE_INTERVAL_TREE] = {
      [BF_TYPE_TREE] =          true,
      [BF_TYPE_INTERVAL_TREE] = true
    },
    [BF_TYPE_OCTREE] = {
      [BF_TYPE_TREE] =   true,
      [BF_TYPE_OCTREE] = true
    },
    [BF_TYPE_QUADTREE] = {
      [BF_TYPE_TREE] =     true,
      [BF_TYPE_QUADTREE] = true
    },

    /** Tree node hierarchy* */
    [BF_TYPE_TREE_NODE] = {
      [BF_TYPE_TREE_NODE] = true
    },
    [BF_TYPE_FIEDLER_TREE_NODE] = {
      [BF_TYPE_TREE_NODE] =         true,
      [BF_TYPE_FIEDLER_TREE_NODE] = true
    },
    [BF_TYPE_INTERVAL_TREE_NODE] = {
      [BF_TYPE_TREE_NODE] =          true,
      [BF_TYPE_INTERVAL_TREE_NODE] = true
    },
    [BF_TYPE_OCTREE_NODE] = {
      [BF_TYPE_TREE_NODE] =   true,
      [BF_TYPE_OCTREE_NODE] = true
    },
    [BF_TYPE_QUADTREE_NODE] = {
      [BF_TYPE_TREE_NODE] =     true,
      [BF_TYPE_QUADTREE_NODE] = true
    },

    /** Tree iterator hierarchy: */
    [BF_TYPE_TREE_ITER] = {
      [BF_TYPE_TREE_ITER] = true
    },
    [BF_TYPE_TREE_ITER_POST_ORDER] = {
      [BF_TYPE_TREE_ITER] =            true,
      [BF_TYPE_TREE_ITER_POST_ORDER] = true
    },
  };
  return _[derived][parent];
}
