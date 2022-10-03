#pragma once

#include <stdbool.h>

typedef struct BfMat BfMat;
typedef struct BfMatBlock BfMatBlock;
typedef struct BfMatBlockCoo BfMatBlockCoo;
typedef struct BfMatBlockDense BfMatBlockDense;
typedef struct BfMatBlockDiag BfMatBlockDiag;
typedef struct BfMatCooComplex BfMatCooComplex;
typedef struct BfMatCooReal BfMatCooReal;
typedef struct BfMatDenseComplex BfMatDenseComplex;
typedef struct BfMatDenseReal BfMatDenseReal;
typedef struct BfMatDiagReal BfMatDiagReal;
typedef struct BfMatGivensComplex BfMatGivensComplex;
typedef struct BfMatProduct BfMatProduct;
typedef struct BfMatSum BfMatSum;

typedef struct BfVec BfVec;
typedef struct BfVecComplex BfVecComplex;
typedef struct BfVecReal BfVecReal;

typedef enum BfTypes {
  /* Mat hierarchy */
  BF_TYPE_MAT,
  BF_TYPE_MAT_BLOCK,
  BF_TYPE_MAT_BLOCK_COO,
  BF_TYPE_MAT_BLOCK_DENSE,
  BF_TYPE_MAT_BLOCK_DIAG,
  BF_TYPE_MAT_COO_COMPLEX,
  BF_TYPE_MAT_COO_REAL,
  BF_TYPE_MAT_DENSE_COMPLEX,
  BF_TYPE_MAT_DENSE_REAL,
  BF_TYPE_MAT_DIAG_REAL,
  BF_TYPE_MAT_GIVENS_COMPLEX,
  BF_TYPE_MAT_PRODUCT,
  BF_TYPE_MAT_SUM,

  /* Vec hierarchy */
  BF_TYPE_VEC,
  BF_TYPE_VEC_COMPLEX,
  BF_TYPE_VEC_REAL,

  /* The total number of types */
  BF_TYPE_COUNT,

  /* Placeholder for missing type information */
  BF_TYPE_UNKNOWN
} BfType;

bool bfTypeDerivedFrom(BfType derived, BfType parent);
