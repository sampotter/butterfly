#include <bf/mat_func.h>

#include <assert.h>
#include <stdlib.h>

#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mem.h>

/** Interface: MatFunc */

static BfMatVtable MAT_VTABLE = {
  .GetNumRows = (__typeof__(&bfMatGetNumRows))bfMatFuncGetNumRows,
  .GetNumCols = (__typeof__(&bfMatGetNumRows))bfMatFuncGetNumCols,
  .Mul = (__typeof__(&bfMatMul))bfMatFuncMul,
};

BfSize bfMatFuncGetNumRows(BfMatFunc const *matFunc) {
  BEGIN_ERROR_HANDLING();

  BfSize numRows = BF_SIZE_BAD_VALUE;

  BfMat const *mat = bfMatFuncConstToMatConst(matFunc);

  if (bfMatIsTransposed(mat))
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);

  numRows = matFunc->super.numRows;

  END_ERROR_HANDLING() {
    assert(false);
  }

  return numRows;
}

BfSize bfMatFuncGetNumCols(BfMatFunc const *matFunc) {
  BEGIN_ERROR_HANDLING();

  BfSize numCols = BF_SIZE_BAD_VALUE;

  BfMat const *mat = bfMatFuncConstToMatConst(matFunc);

  if (bfMatIsTransposed(mat))
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);

  numCols = matFunc->super.numCols;

  END_ERROR_HANDLING() {
    assert(false);
  }

  return numCols;
}

BfMat *bfMatFuncMul(BfMatFunc const *matFunc, BfMat const *otherMat) {
  BEGIN_ERROR_HANDLING();

  if (bfMatFuncGetNumCols(matFunc) != bfMatGetNumRows(otherMat))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfMat *result = matFunc->matMul(otherMat, matFunc->aux);
  HANDLE_ERROR();

  if (bfMatGetNumRows(result) != bfMatFuncGetNumRows(matFunc))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (bfMatGetNumCols(result) != bfMatGetNumCols(otherMat))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  END_ERROR_HANDLING() {
    assert(false);
  }

  return result;
}

/** Upcasting: MatFunc -> Mat */

BfMat *bfMatFuncToMat(BfMatFunc *matFunc) {
  return &matFunc->super;
}

BfMat const *bfMatFuncConstToMatConst(BfMatFunc const *matFunc) {
  return &matFunc->super;
}

/** Implementation: MatFunc */

BfMatFunc *bfMatFuncNew() {
  BEGIN_ERROR_HANDLING();

  BfMatFunc *matFunc = bfMemAlloc(1, sizeof(BfMatFunc));
  if (matFunc == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING() {}

  return matFunc;
}

void bfMatFuncInit(BfMatFunc *matFunc, BfSize numRows, BfSize numCols, MatMulFunc matMul, void *aux) {
  bfMatInit(&matFunc->super, &MAT_VTABLE, numRows, numCols);

  matFunc->matMul = matMul;
  matFunc->aux = aux;
}

void bfMatFuncDeinit(BfMatFunc *matFunc) {}

void bfMatFuncDealloc(BfMatFunc **matFunc) {
  free(*matFunc);
  *matFunc = NULL;
}

void bfMatFuncDeinitAndDealloc(BfMatFunc **matFunc) {
  bfMatFuncDeinit(*matFunc);
  bfMatFuncDealloc(matFunc);
}
