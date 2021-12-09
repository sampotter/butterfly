#include "dtype.h"

#include <assert.h>

#include "mat.h"

BfSize bfDtypeSize(enum BfDtypes dtype) {
  static BfSize size[] = {
    sizeof(BfReal),
    sizeof(BfComplex),
    sizeof(BfMat)
  };
  return size[dtype];
}

enum BfError bfSizeOfDtype(enum BfDtypes dtype, BfSize *nbytes) {
  switch (dtype) {
  case BF_DTYPE_REAL:
    *nbytes = sizeof(BfReal);
    return BF_ERROR_NO_ERROR;
  case BF_DTYPE_COMPLEX:
    *nbytes = sizeof(BfComplex);
    return BF_ERROR_NO_ERROR;
  case BF_DTYPE_MAT:
    *nbytes = sizeof(BfMat);
    return BF_ERROR_NO_ERROR;
  default:
    return BF_ERROR_INVALID_ARGUMENTS;
  };
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
