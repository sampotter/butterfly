#include <bf/mat_types.h>

bool bfMatTypeDerivedFrom(BfMatType derived, BfMatType parent) {
  static bool _[BF_MAT_TYPE_COUNT][BF_MAT_TYPE_COUNT] = {
    [BF_MAT_TYPE_MAT] = {
      [BF_MAT_TYPE_MAT]           = true,
      [BF_MAT_TYPE_BLOCK]         = false,
      [BF_MAT_TYPE_BLOCK_COO]     = false,
      [BF_MAT_TYPE_BLOCK_DENSE]   = false,
      [BF_MAT_TYPE_BLOCK_DIAG]    = false,
      [BF_MAT_TYPE_DENSE_COMPLEX] = false,
      [BF_MAT_TYPE_DIAG_REAL]     = false,
      [BF_MAT_TYPE_PRODUCT]       = false
    },
    [BF_MAT_TYPE_BLOCK] = {
      [BF_MAT_TYPE_MAT]           = true,
      [BF_MAT_TYPE_BLOCK]         = true,
      [BF_MAT_TYPE_BLOCK_COO]     = false,
      [BF_MAT_TYPE_BLOCK_DENSE]   = false,
      [BF_MAT_TYPE_BLOCK_DIAG]    = false,
      [BF_MAT_TYPE_DENSE_COMPLEX] = false,
      [BF_MAT_TYPE_DIAG_REAL]     = false,
      [BF_MAT_TYPE_PRODUCT]       = false
    },
    [BF_MAT_TYPE_BLOCK_COO] = {
      [BF_MAT_TYPE_MAT]           = true,
      [BF_MAT_TYPE_BLOCK]         = true,
      [BF_MAT_TYPE_BLOCK_COO]     = true,
      [BF_MAT_TYPE_BLOCK_DENSE]   = false,
      [BF_MAT_TYPE_BLOCK_DIAG]    = false,
      [BF_MAT_TYPE_DENSE_COMPLEX] = false,
      [BF_MAT_TYPE_DIAG_REAL]     = false,
      [BF_MAT_TYPE_PRODUCT]       = false
    },
    [BF_MAT_TYPE_BLOCK_DENSE] = {
      [BF_MAT_TYPE_MAT]           = true,
      [BF_MAT_TYPE_BLOCK]         = true,
      [BF_MAT_TYPE_BLOCK_COO]     = false,
      [BF_MAT_TYPE_BLOCK_DENSE]   = true,
      [BF_MAT_TYPE_BLOCK_DIAG]    = false,
      [BF_MAT_TYPE_DENSE_COMPLEX] = false,
      [BF_MAT_TYPE_DIAG_REAL]     = false,
      [BF_MAT_TYPE_PRODUCT]       = false
    },
    [BF_MAT_TYPE_BLOCK_DIAG] = {
      [BF_MAT_TYPE_MAT]           = true,
      [BF_MAT_TYPE_BLOCK]         = true,
      [BF_MAT_TYPE_BLOCK_COO]     = false,
      [BF_MAT_TYPE_BLOCK_DENSE]   = false,
      [BF_MAT_TYPE_BLOCK_DIAG]    = true,
      [BF_MAT_TYPE_DENSE_COMPLEX] = false,
      [BF_MAT_TYPE_DIAG_REAL]     = false,
      [BF_MAT_TYPE_PRODUCT]       = false
    },
    [BF_MAT_TYPE_DENSE_COMPLEX] = {
      [BF_MAT_TYPE_MAT]           = true,
      [BF_MAT_TYPE_BLOCK]         = false,
      [BF_MAT_TYPE_BLOCK_COO]     = false,
      [BF_MAT_TYPE_BLOCK_DENSE]   = false,
      [BF_MAT_TYPE_BLOCK_DIAG]    = false,
      [BF_MAT_TYPE_DENSE_COMPLEX] = true,
      [BF_MAT_TYPE_DIAG_REAL]     = false,
      [BF_MAT_TYPE_PRODUCT]       = false
    },
    [BF_MAT_TYPE_DIAG_REAL] = {
      [BF_MAT_TYPE_MAT]           = true,
      [BF_MAT_TYPE_BLOCK]         = false,
      [BF_MAT_TYPE_BLOCK_COO]     = false,
      [BF_MAT_TYPE_BLOCK_DENSE]   = false,
      [BF_MAT_TYPE_BLOCK_DIAG]    = false,
      [BF_MAT_TYPE_DENSE_COMPLEX] = false,
      [BF_MAT_TYPE_DIAG_REAL]     = true,
      [BF_MAT_TYPE_PRODUCT]       = false
    },
    [BF_MAT_TYPE_PRODUCT] = {
      [BF_MAT_TYPE_MAT]           = true,
      [BF_MAT_TYPE_BLOCK]         = false,
      [BF_MAT_TYPE_BLOCK_COO]     = false,
      [BF_MAT_TYPE_BLOCK_DENSE]   = false,
      [BF_MAT_TYPE_BLOCK_DIAG]    = false,
      [BF_MAT_TYPE_DENSE_COMPLEX] = false,
      [BF_MAT_TYPE_DIAG_REAL]     = false,
      [BF_MAT_TYPE_PRODUCT]       = true
    }
  };
  return _[derived][parent];
}
