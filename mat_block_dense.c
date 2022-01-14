#include "mat_block_dense.h"

#include <assert.h>
#include <stdlib.h>

#include "error.h"
#include "error_macros.h"

static BfMatVtable matBlockDenseVtbl = {
  .deinit = (__typeof__(&bfMatDeinit))bfMatBlockDenseDeinit,
  .delete = (__typeof__(&bfMatDelete))bfMatBlockDenseDelete,
  .deinitAndDelete = (__typeof__(&bfMatDeinitAndDelete))bfMatBlockDenseDeinitAndDelete,
  .getType = (__typeof__(&bfMatGetType))bfMatBlockDenseGetType,
  .numBytes = (__typeof__(&bfMatNumBytes))bfMatBlockDenseNumBytes,
  .save = (__typeof__(&bfMatSave))bfMatBlockDenseSave,
  .mul = (__typeof__(&bfMatMul))bfMatBlockDenseMul,
  .lstSq = (__typeof__(&bfMatLstSq))bfMatBlockDenseLstSq
};

BfMatBlockDense *bfMatBlockDenseNew() {
  BEGIN_ERROR_HANDLING();

  BfMatBlockDense *mat = malloc(sizeof(BfMatBlockDense));
  if (mat == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING() {}

  return mat;
}

void BfMatBlockInit(BfMatBlockDense *mat, BfSize numRows, BfSize numCols) {
  BEGIN_ERROR_HANDLING();

  bfMatInit(&mat->super, &matBlockDenseVtbl, numRows, numCols);

  BfSize numElts = numRows*numCols;

  mat->data = malloc(numElts*sizeof(BfMat *));
  if (mat->data == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  for (BfSize i = 0; i < numElts; ++i)
    mat->data[i] = NULL;

  END_ERROR_HANDLING() {}
}

void bfMatBlockDenseDeinit(BfMatBlockDense *mat) {
  free(mat->data);
  mat->data = NULL;
}

void bfMatBlockDenseDelete(BfMatBlockDense **mat) {
  free(*mat);
  *mat = NULL;
}

void bfMatBlockDenseDeinitAndDelete(BfMatBlockDense **mat) {
  bfMatBlockDenseDeinit(*mat);
  bfMatBlockDenseDelete(mat);
}

BfMat *bfMatBlockDenseGetMatPtr(BfMatBlockDense *mat) {
  return &mat->super;
}

BfMatType bfMatBlockDenseGetType(BfMatBlockDense const *mat) {
  (void)mat;
  return BF_MAT_TYPE_BLOCK_DENSE;
}

BfSize bfMatBlockDenseNumBytes(BfMatBlockDense const *mat) {
  (void)mat;
  assert(false);
  return BF_SIZE_BAD_VALUE;
}

void bfMatBlockDenseSave(BfMatBlockDense const *mat, char const *path) {
  (void)mat;
  (void)path;
  assert(false);
}

BfMat *bfMatBlockDenseMul(BfMatBlockDense const *op1, BfMat const *op2) {
  (void)op1;
  (void)op2;
  assert(false);
  return NULL;
}

BfMat *bfMatBlockDenseLstSq(BfMatBlockDense const *lhs, BfMat const *rhs) {
  (void)lhs;
  (void)rhs;
  assert(false);
  return NULL;
}

void bfMatBlockDenseSetBlock(BfMatBlockDense *mat, BfSize i, BfSize j, BfMat *block) {
  BEGIN_ERROR_HANDLING();

  BfSize m = mat->super.numRows;
  BfSize n = mat->super.numCols;

  if (i >= m)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (j >= n)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  mat->data[n*i + j] = block;

  END_ERROR_HANDLING() {}
}
