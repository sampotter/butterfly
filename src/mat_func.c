#include <bf/mat_func.h>

#include <bf/assert.h>
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
  BF_ERROR_BEGIN();

  BfSize numRows = BF_SIZE_BAD_VALUE;

  BfMat const *mat = bfMatFuncConstToMatConst(matFunc);

  if (bfMatIsTransposed(mat))
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);

  numRows = matFunc->super.numRows;

  BF_ERROR_END() {
    BF_DIE();
  }

  return numRows;
}

BfSize bfMatFuncGetNumCols(BfMatFunc const *matFunc) {
  BF_ERROR_BEGIN();

  BfSize numCols = BF_SIZE_BAD_VALUE;

  BfMat const *mat = bfMatFuncConstToMatConst(matFunc);

  if (bfMatIsTransposed(mat))
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);

  numCols = matFunc->super.numCols;

  BF_ERROR_END() {
    BF_DIE();
  }

  return numCols;
}

BfMat *bfMatFuncMul(BfMatFunc const *matFunc, BfMat const *otherMat) {
  BF_ERROR_BEGIN();

  if (bfMatFuncGetNumCols(matFunc) != bfMatGetNumRows(otherMat))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfMat *result = matFunc->matMul(otherMat, matFunc->aux);
  HANDLE_ERROR();

  if (bfMatGetNumRows(result) != bfMatFuncGetNumRows(matFunc))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (bfMatGetNumCols(result) != bfMatGetNumCols(otherMat))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BF_ERROR_END() {
    BF_DIE();
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
  BF_ERROR_BEGIN();

  BfMatFunc *matFunc = bfMemAlloc(1, sizeof(BfMatFunc));
  if (matFunc == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BF_ERROR_END() {}

  return matFunc;
}

void bfMatFuncInit(BfMatFunc *matFunc, BfSize numRows, BfSize numCols, MatMulFunc matMul, void *aux) {
  bfMatInit(&matFunc->super, &MAT_VTABLE, numRows, numCols);

  matFunc->matMul = matMul;
  matFunc->aux = aux;
}

void bfMatFuncDeinit(BfMatFunc *matFunc) {}

void bfMatFuncDealloc(BfMatFunc **matFunc) {
  bfMemFree(*matFunc);
  *matFunc = NULL;
}

void bfMatFuncDeinitAndDealloc(BfMatFunc **matFunc) {
  bfMatFuncDeinit(*matFunc);
  bfMatFuncDealloc(matFunc);
}
