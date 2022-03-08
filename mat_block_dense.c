#include "mat_block_dense.h"

#include <assert.h>
#include <stdlib.h>

#include "error.h"
#include "error_macros.h"

static BfMatVtable matVtbl = {
  .deinit = (__typeof__(&bfMatDeinit))bfMatBlockDenseDeinit,
  .delete = (__typeof__(&bfMatDelete))bfMatBlockDenseDelete,
  .deinitAndDelete = (__typeof__(&bfMatDeinitAndDelete))bfMatBlockDenseDeinitAndDelete,
  .getType = (__typeof__(&bfMatGetType))bfMatBlockDenseGetType,
  .numBytes = (__typeof__(&bfMatNumBytes))bfMatBlockDenseNumBytes,
  .save = (__typeof__(&bfMatSave))bfMatBlockDenseSave,
  .mul = (__typeof__(&bfMatMul))bfMatBlockDenseMul,
  .lstSq = (__typeof__(&bfMatLstSq))bfMatBlockDenseLstSq
};

static BfMatBlockVtable matBlockVtbl = {
  .numBlocks = (__typeof__(&bfMatBlockNumBlocks))bfMatBlockDenseNumBlocks
};

BfMatBlockDense *bfMatBlockDenseNew() {
  return malloc(sizeof(BfMatBlockDense));
}

void bfMatBlockDenseInit(BfMatBlockDense *mat,
                         BfSize numBlockRows, BfSize numBlockCols) {
  BEGIN_ERROR_HANDLING();

  BfSize numBlocks = numBlockRows*numBlockCols;

  bfMatBlockInit(&mat->super,
                 &matVtbl, &matBlockVtbl,
                 numBlocks, numBlockRows, numBlockCols);
  HANDLE_ERROR();

  END_ERROR_HANDLING()
    bfMatBlockDeinit(&mat->super);
}

BfMat *bfMatBlockDenseGetMatPtr(BfMatBlockDense *mat) {
  return &mat->super.super;
}

void bfMatBlockDenseDeinit(BfMatBlockDense *mat) {
  (void)mat;
  assert(false);
}

void bfMatBlockDenseDelete(BfMatBlockDense **mat) {
  free(*mat);
  *mat = NULL;
}

void bfMatBlockDenseDeinitAndDelete(BfMatBlockDense **mat) {
  bfMatBlockDenseDeinit(*mat);
  bfMatBlockDenseDelete(mat);
}

BfMatType bfMatBlockDenseGetType(BfMatBlockDense *mat) {
  (void)mat;
  return BF_MAT_TYPE_BLOCK_DENSE;
}

BfSize bfMatBlockDenseNumBytes(BfMatBlockDense *mat) {
  BfSize num_bytes = 0;

  /* memory occupied by mat itself */
  num_bytes += sizeof(BfMatBlockDense);

  /* memory occupied by blocks */
  BfSize numBlocks = bfMatBlockDenseNumBlocks(mat);
  for (BfSize i = 0; i < numBlocks; ++i)
    num_bytes += bfMatNumBytes(mat->super.block[i]);

  /* mat->super.rowOffset and mat->super.colOffset */
  num_bytes += (mat->super.super.numRows + 1)*sizeof(BfSize);
  num_bytes += (mat->super.super.numCols + 1)*sizeof(BfSize);

  /* mat->rowInd and mat->colInd */
  num_bytes += 2*numBlocks*sizeof(BfSize);

  return num_bytes;
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

  // TODO: not sure whether this should be implemented or not... see:
  //   https://en.wikipedia.org/wiki/Block_matrix_pseudoinverse
}

BfSize bfMatBlockDenseNumBlocks(BfMatBlockDense const *mat) {
  BfSize numBlockRows = mat->super.super.numRows;
  BfSize numBlockCols = mat->super.super.numCols;
  return numBlockRows*numBlockCols;
}
