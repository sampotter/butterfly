#include "mat_block_dense.h"

#include <assert.h>
#include <stdlib.h>

#include "error.h"
#include "error_macros.h"

BF_DEFINE_MAT_VTABLE(MatBlockDense);

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

BfMat *bfMatBlockDenseGetBlock(BfMatBlockDense *mat, BfSize i, BfSize j) {
  BEGIN_ERROR_HANDLING();

  BfMat *block = NULL;

  BfSize numBlockRows = mat->super.super.numRows;
  if (i >= numBlockRows)
    RAISE_ERROR(BF_ERROR_OUT_OF_RANGE);

  BfSize numBlockCols = mat->super.super.numCols;
  if (j >= numBlockCols)
    RAISE_ERROR(BF_ERROR_OUT_OF_RANGE);

  block = mat->super.block[numBlockCols*i + j];

  END_ERROR_HANDLING() {}

  return block;
}

void bfMatBlockDenseSetBlock(BfMatBlockDense *mat, BfSize i, BfSize j,
                             BfMat *block) {
  BEGIN_ERROR_HANDLING();

  BfSize numBlockRows = mat->super.super.numRows;
  if (i >= numBlockRows)
    RAISE_ERROR(BF_ERROR_OUT_OF_RANGE);

  BfSize numBlockCols = mat->super.super.numCols;
  if (j >= numBlockCols)
    RAISE_ERROR(BF_ERROR_OUT_OF_RANGE);

  mat->super.block[numBlockCols*i + j] = block;

  END_ERROR_HANDLING() {}
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

BfMat *bfMatBlockDenseZerosLike(BfMatBlockDense *mat, BfSize numRows, BfSize numCols) {
  (void)mat;
  assert(false);
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

BfSize bfMatBlockDenseGetNumRows(BfMatBlockDense const *mat) {
  return mat->super.rowOffset[mat->super.super.numRows];
}

BfSize bfMatBlockDenseGetNumCols(BfMatBlockDense const *mat) {
  return mat->super.colOffset[mat->super.super.numCols];
}

BfMatBlockDense *bfMatBlockDenseGetRowRange(BfMatBlockDense *mat, BfSize i0, BfSize i1) {
  BEGIN_ERROR_HANDLING();

  BfMatBlockDense *view = NULL;

  if (i0 >= i1)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize numRows = bfMatBlockDenseGetNumRows(mat);

  if (i0 >= numRows)
    RAISE_ERROR(BF_ERROR_OUT_OF_RANGE);

  if (i1 >= numRows)
    RAISE_ERROR(BF_ERROR_OUT_OF_RANGE);

  /* Find the offset of each row of blocks. For now, we require that
   * i0 and i1 exactly correspond to the start of two distinct block
   * rows.
   *
   * In the future, we could relax this assumption and allow block
   * matrices to get partitioned when i0 and i1 don't exactly align
   * with the block partitions. But probably it would be better to
   * just do that in a separate function. */

  BfSize p0 = bfMatBlockFindRowBlock(&mat->super, i0);
  if (p0 == BF_SIZE_BAD_VALUE || i0 != mat->super.rowOffset[p0])
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize p1 = bfMatBlockFindRowBlock(&mat->super, i1);
  if (p1 == BF_SIZE_BAD_VALUE || i1 != mat->super.rowOffset[p1])
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize numColBlocks = bfMatBlockGetNumColBlocks(&mat->super);

  view = bfMatBlockDenseNew();
  HANDLE_ERROR();

  bfMatBlockDenseInit(view, p1 - p0, numColBlocks);

  /* the row range is a view of the original matrix */
  bfMatBlockDenseGetMatPtr(view)->props |= BF_MAT_PROPS_VIEW;

  /* point to the subrange of blocks */
  view->super.block = mat->super.block + numColBlocks*p0;

  END_ERROR_HANDLING() {}

  return view;
}

BfMatBlockDense *bfMatBlockDenseGetColRange(BfMatBlockDense *mat, BfSize j0, BfSize j1) {
  (void)mat;
  (void)j0;
  (void)j1;
  assert(false);
  return NULL;
}

void *bfMatBlockDenseAddInplace(BfMatBlockDense *op1, BfMat *op2) {
  (void)op1;
  (void)op2;
  assert(false);
}

BfMat *bfMatBlockDenseMul(BfMatBlockDense const *op1, BfMat const *op2) {
  BEGIN_ERROR_HANDLING();

  BfSize numRows, numCols, numRowBlocks, numColBlocks;
  BfMatBlock const *super1 = &op1->super;
  BfMat *block = NULL, *op2Rows = NULL;
  BfMat *result = NULL, *resultRows = NULL, *tmp = NULL;

  numRowBlocks = bfMatBlockGetNumRowBlocks(super1);
  numColBlocks = bfMatBlockGetNumColBlocks(super1);

  numRows = bfMatGetNumRows(&super1->super);
  numCols = bfMatGetNumCols(op2);

  result = bfMatZerosLike(op2, numRows, numCols);
  if (result == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  for (BfSize i = 0, i0, i1; i < numRowBlocks; ++i) {
    i0 = super1->rowOffset[i];
    i1 = super1->rowOffset[i + 1];
    resultRows = bfMatGetRowRange(result, i0, i1);
    for (BfSize j = 0, j0, j1; j < numColBlocks; ++j) {
      j0 = super1->colOffset[j];
      j1 = super1->colOffset[j + 1];
      op2Rows = bfMatGetRowRange((BfMat *)op2, j0, j1);
      block = bfMatBlockDenseGetBlock((BfMatBlockDense *)op1, i, j);
      tmp = bfMatMul(block, op2Rows);
      bfMatAddInplace(resultRows, tmp);
      bfMatDeinitAndDelete(&tmp);
    }
  }

  END_ERROR_HANDLING()
    bfMatDeinitAndDelete(&result);

  return result;
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
