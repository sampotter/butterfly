#include <bf/mat_types.h>

bool bfMatTypeDerivedFrom(BfMatType derived, BfMatType parent) {
  /* For each entry in this lookup table, we should have:
   *
   *   `_[SUBTYPE][TYPE] == true` if `SUBTYPE` is a subtype of `TYPE`
   *
   * When a new type is added, this lookup table needs to be updated.
   *
   * NOTE: static arrays are zero-initialized by default in C99. */
  static bool _[BF_MAT_TYPE_COUNT][BF_MAT_TYPE_COUNT] = {
    [BF_MAT_TYPE_MAT] = {
      [BF_MAT_TYPE_MAT]           = true,
    },
    [BF_MAT_TYPE_BLOCK] = {
      [BF_MAT_TYPE_MAT]           = true,
      [BF_MAT_TYPE_BLOCK]         = true,
    },
    [BF_MAT_TYPE_BLOCK_COO] = {
      [BF_MAT_TYPE_MAT]           = true,
      [BF_MAT_TYPE_BLOCK]         = true,
      [BF_MAT_TYPE_BLOCK_COO]     = true,
    },
    [BF_MAT_TYPE_BLOCK_DENSE] = {
      [BF_MAT_TYPE_MAT]           = true,
      [BF_MAT_TYPE_BLOCK]         = true,
      [BF_MAT_TYPE_BLOCK_DENSE]   = true,
    },
    [BF_MAT_TYPE_BLOCK_DIAG] = {
      [BF_MAT_TYPE_MAT]           = true,
      [BF_MAT_TYPE_BLOCK]         = true,
      [BF_MAT_TYPE_BLOCK_DIAG]    = true,
    },
    [BF_MAT_TYPE_DENSE_COMPLEX] = {
      [BF_MAT_TYPE_MAT]           = true,
      [BF_MAT_TYPE_DENSE_COMPLEX] = true,
    },
    [BF_MAT_TYPE_DENSE_REAL] = {
      [BF_MAT_TYPE_MAT]           = true,
      [BF_MAT_TYPE_DENSE_REAL]    = true,
    },
    [BF_MAT_TYPE_DIAG_REAL] = {
      [BF_MAT_TYPE_MAT]           = true,
      [BF_MAT_TYPE_DIAG_REAL]     = true,
    },
    [BF_MAT_TYPE_PRODUCT] = {
      [BF_MAT_TYPE_MAT]           = true,
      [BF_MAT_TYPE_PRODUCT]       = true
    }
  };
  return _[derived][parent];
}
