#include <bf/mat_block_dense.h>

#include <assert.h>
#include <stdlib.h>

#include <bf/error.h>
#include <bf/error_macros.h>

#define INTERFACE BF_INTERFACE_Mat
BF_DEFINE_VTABLE(Mat, MatBlockDense);
#undef INTERFACE

#define INTERFACE BF_INTERFACE_MatBlock
BF_DEFINE_VTABLE(MatBlock, MatBlockDense);
#undef INTERFACE

BfMat *bfMatBlockDenseToMat(BfMatBlockDense *matBlockDense) {
  return &matBlockDense->super.super;
}

BfMat const *bfMatBlockDenseConstToMatConst(BfMatBlockDense const *matBlockDense) {
  return &matBlockDense->super.super;
}

BfMatBlockDense const *bfMatConstToMatBlockDenseConst(BfMat const *mat) {
  if (!bfMatInstanceOf(mat, BF_MAT_TYPE_BLOCK_DENSE)) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return NULL;
  } else {
    return (BfMatBlockDense const *)mat;
  }
}

BfMatBlockDense *bfMatBlockDenseNew() {
  return malloc(sizeof(BfMatBlockDense));
}

void bfMatBlockDenseInit(BfMatBlockDense *mat,
                         BfSize numBlockRows, BfSize numBlockCols) {
  BEGIN_ERROR_HANDLING();

  BfSize numBlocks = numBlockRows*numBlockCols;

  bfMatBlockInit(&mat->super,
                 &MatVtbl, &MatBlockVtbl,
                 numBlocks, numBlockRows, numBlockCols);
  HANDLE_ERROR();

  END_ERROR_HANDLING()
    bfMatBlockDeinit(&mat->super);
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
  bfMatBlockDeinit(&mat->super);
}

void bfMatBlockDenseDealloc(BfMatBlockDense **mat) {
  free(*mat);
  *mat = NULL;
}

void bfMatBlockDenseDeinitAndDealloc(BfMatBlockDense **mat) {
  bfMatBlockDenseDeinit(*mat);
  bfMatBlockDenseDealloc(mat);
}

void bfMatBlockDenseDelete(BfMat **mat) {
  bfMatBlockDenseDeinitAndDealloc((BfMatBlockDense **)mat);
}

BfMat *bfMatBlockDenseEmptyLike(BfMat const *, BfSize, BfSize) {
  assert(false);
}

BfMat *bfMatBlockDenseZerosLike(BfMat const *, BfSize, BfSize) {
  assert(false);
}

BfMatType bfMatBlockDenseGetType(BfMat const *mat) {
  (void)mat;
  return BF_MAT_TYPE_BLOCK_DENSE;
}

bool bfMatBlockDenseInstanceOf(BfMat const *mat, BfMatType matType) {
  return bfMatTypeDerivedFrom(bfMatGetType(mat), matType);
}

BfSize bfMatBlockDenseNumBytes(BfMat const *mat) {
  BfSize num_bytes = 0;

  BfMatBlock const *matBlock = bfMatConstToMatBlockConst(mat);

  /* memory occupied by mat itself */
  num_bytes += sizeof(BfMatBlockDense);

  /* memory occupied by blocks */
  BfSize numBlocks = bfMatBlockDenseNumBlocks(matBlock);
  for (BfSize i = 0; i < numBlocks; ++i)
    num_bytes += bfMatNumBytes(matBlock->block[i]);

  /* mat->super.rowOffset and mat->super.colOffset */
  num_bytes += (mat->numRows + 1)*sizeof(BfSize);
  num_bytes += (mat->numCols + 1)*sizeof(BfSize);

  /* mat->rowInd and mat->colInd */
  num_bytes += 2*numBlocks*sizeof(BfSize);

  return num_bytes;
}

void bfMatBlockDenseSave(BfMat const *mat, char const *path) {
  (void)mat;
  (void)path;
  assert(false);
}

BfSize bfMatBlockDenseGetNumRows(BfMat const *mat) {
  BfMatBlock const *matBlock = bfMatConstToMatBlockConst(mat);
  return bfMatIsTransposed(mat) ?
    matBlock->colOffset[mat->numCols] :
    matBlock->rowOffset[mat->numRows];
}

BfSize bfMatBlockDenseGetNumCols(BfMat const *mat) {
  BfMatBlock const *matBlock = bfMatConstToMatBlockConst(mat);
  return bfMatIsTransposed(mat) ?
    matBlock->rowOffset[mat->numRows] :
    matBlock->colOffset[mat->numCols];
}

BfMat *bfMatBlockDenseGetRowRange(BfMat *mat, BfSize i0, BfSize i1) {
  BEGIN_ERROR_HANDLING();

  BfMatBlock *matBlock = bfMatToMatBlock(mat);

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

  BfSize p0 = bfMatBlockFindRowBlock(matBlock, i0);
  if (p0 == BF_SIZE_BAD_VALUE || i0 != matBlock->rowOffset[p0])
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize p1 = bfMatBlockFindRowBlock(matBlock, i1);
  if (p1 == BF_SIZE_BAD_VALUE || i1 != matBlock->rowOffset[p1])
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize numColBlocks = bfMatBlockGetNumColBlocks(matBlock);

  view = bfMatBlockDenseNew();
  HANDLE_ERROR();

  bfMatBlockDenseInit(view, p1 - p0, numColBlocks);

  /* the row range is a view of the original matrix */
  bfMatBlockDenseToMat(view)->props |= BF_MAT_PROPS_VIEW;

  /* point to the subrange of blocks */
  view->super.block = matBlock->block + numColBlocks*p0;

  END_ERROR_HANDLING() {}

  return bfMatBlockDenseToMat(view);
}

BfMat *bfMatBlockDenseGetColRange(BfMat *mat, BfSize j0, BfSize j1) {
  assert(false);
  return NULL;
}

void bfMatBlockDenseSetRowRange(BfMat *, BfSize, BfSize, BfMat const *) {
  assert(false);
}

void bfMatBlockDenseAddInplace(BfMat *, BfMat const *) {
  assert(false);
}

BfMat *bfMatBlockDenseMul(BfMat const *mat, BfMat const *otherMat) {
  BEGIN_ERROR_HANDLING();

  BfMatBlock const *matBlock = bfMatConstToMatBlockConst(mat);
  BfMatBlockDense const *matBlockDense = bfMatConstToMatBlockDenseConst(mat);

  BfSize numRows, numCols, numRowBlocks, numColBlocks;
  BfMat *block = NULL, *op2Rows = NULL;
  BfMat *result = NULL, *resultRows = NULL, *tmp = NULL;

  numRowBlocks = bfMatBlockGetNumRowBlocks(matBlock);
  numColBlocks = bfMatBlockGetNumColBlocks(matBlock);

  numRows = bfMatGetNumRows(mat);
  numCols = bfMatGetNumCols(otherMat);

  result = bfMatZerosLike(otherMat, numRows, numCols);
  if (result == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  for (BfSize i = 0, i0, i1; i < numRowBlocks; ++i) {
    i0 = matBlock->rowOffset[i];
    i1 = matBlock->rowOffset[i + 1];
    resultRows = bfMatGetRowRange(result, i0, i1);
    for (BfSize j = 0, j0, j1; j < numColBlocks; ++j) {
      j0 = matBlock->colOffset[j];
      j1 = matBlock->colOffset[j + 1];
      op2Rows = bfMatGetRowRange((BfMat *)otherMat, j0, j1);
      block = bfMatBlockDenseGetBlock((BfMatBlockDense *)matBlockDense, i, j);
      tmp = bfMatMul(block, op2Rows);
      bfMatAddInplace(resultRows, tmp);
      bfMatDelete(&tmp);
    }
  }

  END_ERROR_HANDLING()
    bfMatDelete(&result);

  return result;
}

BfMat *bfMatBlockDenseLstSq(BfMat const *lhs, BfMat const *rhs) {
  (void)lhs;
  (void)rhs;
  assert(false);
  return NULL;

  // TODO: not sure whether this should be implemented or not... see:
  //   https://en.wikipedia.org/wiki/Block_matrix_pseudoinverse
}

BfSize bfMatBlockDenseNumBlocks(BfMatBlock const *matBlock) {
  BfMat const *mat = bfMatBlockConstToMatConst(matBlock);
  return mat->numRows*mat->numCols;
}

BfSize bfMatBlockDenseGetNumRowBlocks(BfMatBlock const *) { assert(false); }

BfSize bfMatBlockDenseGetNumColBlocks(BfMatBlock const *) { assert(false); }

BfSize bfMatBlockDenseGetNumBlockRows(BfMatBlock const *, BfSize) {
  assert(false);
}

BfSize bfMatBlockDenseGetNumBlockCols(BfMatBlock const *, BfSize) {
  assert(false);
}
