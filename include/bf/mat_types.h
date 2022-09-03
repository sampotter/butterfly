#pragma once

#include <stdbool.h>

// TODO: rename: BF_MAT_TYPE_* -> BF_TYPE_MAT_* (so that first type is
// BF_TYPE_MAT, which is consistent with the rest of the typenames)
typedef enum BfMatTypes {
  BF_MAT_TYPE_MAT,
  BF_MAT_TYPE_BLOCK,
  BF_MAT_TYPE_BLOCK_COO,
  BF_MAT_TYPE_BLOCK_DENSE,
  BF_MAT_TYPE_BLOCK_DIAG,
  BF_MAT_TYPE_DENSE_COMPLEX,
  BF_MAT_TYPE_DIAG_REAL,
  BF_MAT_TYPE_PRODUCT,
  BF_MAT_TYPE_COUNT
} BfMatType;

bool bfMatTypeDerivedFrom(BfMatType derived, BfMatType parent);
