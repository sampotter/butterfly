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
    [BF_TYPE_MAT_COO_COMPLEX] = {
      [BF_TYPE_MAT] =             true,
      [BF_TYPE_MAT_COO_COMPLEX] = true
    },
    [BF_TYPE_MAT_COO_REAL] = {
      [BF_TYPE_MAT] =          true,
      [BF_TYPE_MAT_COO_REAL] = true
    },
    [BF_TYPE_MAT_DENSE_COMPLEX] = {
      [BF_TYPE_MAT]               = true,
      [BF_TYPE_MAT_DENSE_COMPLEX] = true,
    },
    [BF_TYPE_MAT_DENSE_REAL] = {
      [BF_TYPE_MAT]            = true,
      [BF_TYPE_MAT_DENSE_REAL] = true,
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
  };
  return _[derived][parent];
}
