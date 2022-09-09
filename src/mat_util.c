#include <bf/mat_util.h>

#include <bf/error_macros.h>
#include <bf/mat_diag_real.h>

BfMat *bfEye(BfSize numRows, BfSize numCols, BfDtype dtype) {
  switch (dtype) {
  case BF_DTYPE_REAL:
    return bfMatDiagRealToMat(bfMatDiagRealEye(numRows, numCols));
  default:
    bfSetError(BF_ERROR_NOT_IMPLEMENTED);
    return NULL;
  }
}
