#include <bf/mat_zero.h>

#include <stdlib.h>

#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mem.h>
#include <bf/vec_zero.h>

/** Interface: Mat */

static BfMatVtable MAT_VTABLE = {
  .GetRowCopy = (__typeof__(&bfMatZeroGetRowCopy))bfMatZeroGetRowCopy,
  .GetType = (__typeof__(&bfMatZeroGetType))bfMatZeroGetType,
  .GetNumRows = (__typeof__(&bfMatZeroGetNumRows))bfMatZeroGetNumRows,
  .GetNumCols = (__typeof__(&bfMatZeroGetNumCols))bfMatZeroGetNumCols,
};

BfVec *bfMatZeroGetRowCopy(BfMat const *mat, BfSize i) {
  BEGIN_ERROR_HANDLING();

  BfSize numRows = bfMatGetNumRows(mat);
  if (i >= numRows)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize numCols = bfMatGetNumCols(mat);

  if (!bfMatInstanceOf(mat, BF_TYPE_MAT_ZERO))
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  BfVecZero *rowCopy = bfVecZeroNew();
  HANDLE_ERROR();

  bfVecZeroInit(rowCopy, numCols);

  END_ERROR_HANDLING() {}

  return bfVecZeroToVec(rowCopy);
}

BfType bfMatZeroGetType(BfMat const *mat) {
  (void)mat;
  return BF_TYPE_MAT_ZERO;
}

BfSize bfMatZeroGetNumRows(BfMat const *mat) {
  BEGIN_ERROR_HANDLING();

  BfSize numRows = BF_SIZE_BAD_VALUE;

  BfMatZero const *matZero = bfMatConstToMatZeroConst(mat);
  HANDLE_ERROR();

  numRows = matZero->super.numRows;

  END_ERROR_HANDLING() {}

  return numRows;
}

BfSize bfMatZeroGetNumCols(BfMat const *mat) {
  BEGIN_ERROR_HANDLING();

  BfSize numCols = BF_SIZE_BAD_VALUE;

  BfMatZero const *matZero = bfMatConstToMatZeroConst(mat);
  HANDLE_ERROR();

  numCols = matZero->super.numCols;

  END_ERROR_HANDLING() {}

  return numCols;
}

/** Upcasting: */

BfMat *bfMatZeroToMat(BfMatZero *mat) {
  return &mat->super;
}

BfMat const *bfMatZeroConstToMatConst(BfMatZero const *mat) {
  return &mat->super;
}

/** Downcasting: */

BfMatZero const *bfMatConstToMatZeroConst(BfMat const *mat) {
  if (!bfMatInstanceOf(mat, BF_TYPE_MAT_ZERO)) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return NULL;
  } else {
    return (BfMatZero const *)mat;
  }
}

/** Implementation: MatZero */

BfMatZero *bfMatZeroNew() {
  BEGIN_ERROR_HANDLING();

  BfMatZero *mat = bfMemAlloc(1, sizeof(BfMatZero));
  if (mat == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING() {}

  return mat;
}

void bfMatZeroInit(BfMatZero *mat, BfSize numRows, BfSize numCols) {
  BEGIN_ERROR_HANDLING();

  bfMatInit(&mat->super, &MAT_VTABLE, numRows, numCols);
  HANDLE_ERROR();

  END_ERROR_HANDLING()
    bfMatDeinit(&mat->super);
}

void bfMatZeroDeinit(BfMatZero *mat) {
  (void)mat;
}

void bfMatZeroDealloc(BfMatZero **mat) {
  free(*mat);
  *mat = NULL;
}

void bfMatZeroDeinitAndDealloc(BfMatZero **mat) {
  bfMatZeroDeinit(*mat);
  bfMatZeroDealloc(mat);
}
