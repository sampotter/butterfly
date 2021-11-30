#include "dtype.h"

enum BfError bfSizeOfDtype(enum BfDtypes dtype, BfSize *nbytes) {
  switch (dtype) {
  case BF_DTYPE_REAL:
    *nbytes = sizeof(BfReal);
    return BF_ERROR_NO_ERROR;
  case BF_DTYPE_COMPLEX:
    *nbytes = sizeof(BfComplex);
    return BF_ERROR_NO_ERROR;
  default:
    return BF_ERROR_INVALID_ARGUMENTS;
  };
}

bool bfDtypeIsValid(enum BfDtypes dtype) {
  return dtype & (BF_DTYPE_REAL | BF_DTYPE_COMPLEX);
}
