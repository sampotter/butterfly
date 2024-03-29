#include <bf/mat_zero.h>

#include <stdlib.h>

#include <bf/assert.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mem.h>
#include <bf/vec_zero.h>

/** Interface: Mat */

static BfMatVtable MAT_VTABLE = {
  .GetRowCopy = (__typeof__(&bfMatZeroGetRowCopy))bfMatZeroGetRowCopy,
  .Delete = (__typeof__(&bfMatDelete))bfMatZeroDelete,
  .GetType = (__typeof__(&bfMatZeroGetType))bfMatZeroGetType,
  .GetNumRows = (__typeof__(&bfMatZeroGetNumRows))bfMatZeroGetNumRows,
  .GetNumCols = (__typeof__(&bfMatZeroGetNumCols))bfMatZeroGetNumCols,
};

BfVec *bfMatZeroGetRowCopy(BfMat const *mat, BfSize i) {
  BF_ERROR_BEGIN();

  BfVecZero *rowCopy = NULL;

  BfSize numRows = bfMatGetNumRows(mat);
  if (i >= numRows)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize numCols = bfMatGetNumCols(mat);

  if (!bfMatInstanceOf(mat, BF_TYPE_MAT_ZERO))
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  rowCopy = bfVecZeroNew();
  HANDLE_ERROR();

  bfVecZeroInit(rowCopy, numCols);

  BF_ERROR_END() {}

  return bfVecZeroToVec(rowCopy);
}

void bfMatZeroDelete(BfMatZero **matZero) {
  bfMatZeroDeinitAndDealloc(matZero);
}

BfType bfMatZeroGetType(BfMat const *mat) {
  (void)mat;
  return BF_TYPE_MAT_ZERO;
}

BfSize bfMatZeroGetNumRows(BfMat const *mat) {
  BF_ERROR_BEGIN();

  BfSize numRows = BF_SIZE_BAD_VALUE;

  BfMatZero const *matZero = bfMatConstToMatZeroConst(mat);
  HANDLE_ERROR();

  numRows = matZero->super.numRows;

  BF_ERROR_END() {}

  return numRows;
}

BfSize bfMatZeroGetNumCols(BfMat const *mat) {
  BF_ERROR_BEGIN();

  BfSize numCols = BF_SIZE_BAD_VALUE;

  BfMatZero const *matZero = bfMatConstToMatZeroConst(mat);
  HANDLE_ERROR();

  numCols = matZero->super.numCols;

  BF_ERROR_END() {}

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
  BF_ERROR_BEGIN();

  BfMatZero *mat = bfMemAlloc(1, sizeof(BfMatZero));
  if (mat == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BF_ERROR_END() {
    BF_DIE();
  }

  return mat;
}

void bfMatZeroInit(BfMatZero *mat, BfSize numRows, BfSize numCols) {
  BF_ERROR_BEGIN();

  bfMatInit(&mat->super, &MAT_VTABLE, numRows, numCols);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }
}

void bfMatZeroDeinit(BfMatZero *mat) {
  bfMatDeinit(&mat->super);
}

void bfMatZeroDealloc(BfMatZero **mat) {
  free(*mat);
  *mat = NULL;
}

void bfMatZeroDeinitAndDealloc(BfMatZero **mat) {
  bfMatZeroDeinit(*mat);
  bfMatZeroDealloc(mat);
}
