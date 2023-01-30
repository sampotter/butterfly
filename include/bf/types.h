#pragma once

#include <stdbool.h>

typedef struct BfMat BfMat;
typedef struct BfMatBlock BfMatBlock;
typedef struct BfMatBlockCoo BfMatBlockCoo;
typedef struct BfMatBlockDense BfMatBlockDense;
typedef struct BfMatBlockDiag BfMatBlockDiag;
typedef struct BfMatCooComplex BfMatCooComplex;
typedef struct BfMatCooReal BfMatCooReal;
typedef struct BfMatCsrReal BfMatCsrReal;
typedef struct BfMatDense BfMatDense;
typedef struct BfMatDenseComplex BfMatDenseComplex;
typedef struct BfMatDenseReal BfMatDenseReal;
typedef struct BfMatDiagReal BfMatDiagReal;
typedef struct BfMatGivensComplex BfMatGivensComplex;
typedef struct BfMatProduct BfMatProduct;
typedef struct BfMatSum BfMatSum;
typedef struct BfMatZero BfMatZero;

typedef struct BfVec BfVec;
typedef struct BfVecComplex BfVecComplex;
typedef struct BfVecReal BfVecReal;
typedef struct BfVecZero BfVecZero;

typedef struct BfLu BfLu;

typedef struct BfTree BfTree;
typedef struct BfIntervalTree BfIntervalTree;
typedef struct BfQuadtree BfQuadtree;
typedef struct BfOctree BfOctree;
typedef struct BfFiedlerTree BfFiedlerTree;

typedef struct BfTreeNode BfTreeNode;
typedef struct BfIntervalTreeNode BfIntervalTreeNode;
typedef struct BfQuadtreeNode BfQuadtreeNode;
typedef struct BfOctreeNode BfOctreeNode;
typedef struct BfFiedlerTreeNode BfFiedlerTreeNode;

typedef struct BfTreeIter BfTreeIter;
typedef struct BfTreeIterPostOrder BfTreeIterPostOrder;

typedef enum BfTypes {
  /* Mat hierarchy */
  BF_TYPE_MAT,
  BF_TYPE_MAT_COO_COMPLEX,
  BF_TYPE_MAT_COO_REAL,
  BF_TYPE_MAT_CSR_REAL,
  BF_TYPE_MAT_DIAG_REAL,
  BF_TYPE_MAT_GIVENS_COMPLEX,
  BF_TYPE_MAT_PRODUCT,
  BF_TYPE_MAT_SUM,
  BF_TYPE_MAT_ZERO,

  /* MatBlock hierarchy: */
  BF_TYPE_MAT_BLOCK,
  BF_TYPE_MAT_BLOCK_COO,
  BF_TYPE_MAT_BLOCK_DENSE,
  BF_TYPE_MAT_BLOCK_DIAG,

  /* MatDense hierarchy: */
  BF_TYPE_MAT_DENSE,
  BF_TYPE_MAT_DENSE_COMPLEX,
  BF_TYPE_MAT_DENSE_REAL,

  /* Vec hierarchy */
  BF_TYPE_VEC,
  BF_TYPE_VEC_COMPLEX,
  BF_TYPE_VEC_REAL,
  BF_TYPE_VEC_ZERO,

  /* Tree hierarchy */
  BF_TYPE_TREE,
  BF_TYPE_FIEDLER_TREE,
  BF_TYPE_INTERVAL_TREE,
  BF_TYPE_OCTREE,
  BF_TYPE_QUADTREE,

  /* Tree node hierarchy */
  BF_TYPE_TREE_NODE,
  BF_TYPE_FIEDLER_TREE_NODE,
  BF_TYPE_INTERVAL_TREE_NODE,
  BF_TYPE_OCTREE_NODE,
  BF_TYPE_QUADTREE_NODE,

  /* Tree iterator hierarchy */
  BF_TYPE_TREE_ITER,
  BF_TYPE_TREE_ITER_POST_ORDER,

  /* The total number of types */
  BF_TYPE_COUNT,

  /* Placeholder for missing type information */
  BF_TYPE_UNKNOWN
} BfType;

bool bfTypeDerivedFrom(BfType derived, BfType parent);
