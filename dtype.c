#include "dtype.h"

#include <assert.h>

#include "mat.h"

BfSize bfDtypeSize(enum BfDtypes dtype) {
  static BfSize size[] = {
    [BF_DTYPE_VOID]    = 0,
    [BF_DTYPE_REAL]    = sizeof(BfReal),
    [BF_DTYPE_COMPLEX] = sizeof(BfComplex),
    [BF_DTYPE_MAT]     = sizeof(BfMat)
  };
  return size[dtype];
}

static enum BfDtypes const DTYPE_VALID_MASK =
  BF_DTYPE_REAL |
  BF_DTYPE_COMPLEX |
  BF_DTYPE_MAT;

bool bfDtypeIsValid(enum BfDtypes dtype) {
  return dtype & DTYPE_VALID_MASK;
}

void bfGetDtypeZero(enum BfDtypes dtype, BfPtr ptr) {
  assert(bfDtypeIsValid(dtype) && dtype != BF_DTYPE_MAT);

  if (dtype == BF_DTYPE_REAL)    *(BfReal *)ptr = 0;
  if (dtype == BF_DTYPE_COMPLEX) *(BfComplex *)ptr = 0;
}
