#include <bf/mat_block_coo.h>

#include <assert.h>
#include <stdlib.h>

#include <bf/error.h>
#include <bf/error_macros.h>

BF_DEFINE_MAT_VTABLE(MatBlockCoo);

static BfMatBlockVtable matBlockVtbl = {
  .numBlocks = (__typeof__(&bfMatBlockNumBlocks))bfMatBlockCooNumBlocks
};

BfMatBlockCoo *bfMatBlockCooNew() {
  return malloc(sizeof(BfMatBlockCoo));
}

void bfMatBlockCooInit(BfMatBlockCoo *mat, BfSize numBlockRows,
                       BfSize numBlockCols, BfSize numBlocks)
{
  BEGIN_ERROR_HANDLING();

  bfMatBlockInit(&mat->super,
                 &matVtbl, &matBlockVtbl,
                 numBlocks, numBlockRows, numBlockCols);
  HANDLE_ERROR();

  mat->numBlocks = numBlocks;

  /* allocate and initialize row indices to dummy values */
  mat->rowInd = malloc(numBlocks*sizeof(BfSize));
  if (mat->rowInd == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);
  for (BfSize i = 0; i < numBlocks; ++i)
    mat->rowInd[i] = BF_SIZE_BAD_VALUE;

  /* allocate and initialize column indices to dummy values */
  mat->colInd = malloc(numBlocks*sizeof(BfSize));
  if (mat->colInd == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);
  for (BfSize i = 0; i < numBlocks; ++i)
    mat->colInd[i] = BF_SIZE_BAD_VALUE;

  END_ERROR_HANDLING() {
    bfMatBlockDeinit(&mat->super);
    free(mat->rowInd);
    free(mat->colInd);
  }
}

BfMat *bfMatBlockCooGetMatPtr(BfMatBlockCoo *mat) {
  return &mat->super.super;
}

BfMat const *bfMatBlockCooGetMatConstPtr(BfMatBlockCoo const *mat) {
  return &mat->super.super;
}

BfMatBlockCoo *bfMatBlockCooEmptyLike(BfMatBlockCoo const *, BfSize, BfSize) {
  assert(false);
  return NULL;
}

BfMatBlockCoo *bfMatBlockCooZerosLike(BfMatBlockCoo const *, BfSize, BfSize) {
  assert(false);
  return NULL;
}

void bfMatBlockCooDeinit(BfMatBlockCoo *mat) {
  (void)mat;
  assert(false);
}

void bfMatBlockCooDelete(BfMatBlockCoo **mat) {
  free(*mat);
  *mat = NULL;
}

void bfMatBlockCooDeinitAndDelete(BfMatBlockCoo **mat) {
  bfMatBlockCooDeinit(*mat);
  bfMatBlockCooDelete(mat);
}

BfMatType bfMatBlockCooGetType(BfMatBlockCoo const *mat) {
  (void)mat;
  return BF_MAT_TYPE_BLOCK_COO;
}

bool bfMatBlockCooInstanceOf(BfMatBlockCoo const *mat, BfMatType matType) {
  BfMat const *parent = bfMatBlockCooGetMatConstPtr(mat);
  return bfMatTypeDerivedFrom(bfMatGetType(parent), matType);
}

BfSize bfMatBlockCooNumBytes(BfMatBlockCoo const *mat) {
  BfSize num_bytes = 0;

  /* memory occupied by mat itself */
  num_bytes += sizeof(BfMatBlockCoo);

  /* memory occupied by blocks */
  for (BfSize i = 0; i < mat->numBlocks; ++i)
    num_bytes += bfMatNumBytes(mat->super.block[i]);

  /* mat->super.rowOffset and mat->super.colOffset */
  num_bytes += (mat->super.super.numRows + 1)*sizeof(BfSize);
  num_bytes += (mat->super.super.numCols + 1)*sizeof(BfSize);

  /* mat->rowInd and mat->colInd */
  num_bytes += 2*mat->numBlocks*sizeof(BfSize);

  return num_bytes;
}

void bfMatBlockCooSave(BfMatBlockCoo const *, char const *) {
  assert(false);
}

BfSize bfMatBlockCooGetNumRows(BfMatBlockCoo const *mat) {
  BfMatBlock const *super = &mat->super;
  return bfMatIsTransposed(&super->super) ?
    super->colOffset[bfMatBlockGetNumColBlocks(super)] :
    super->rowOffset[bfMatBlockGetNumRowBlocks(super)];
}

BfSize bfMatBlockCooGetNumCols(BfMatBlockCoo const *mat) {
  BfMatBlock const *super = &mat->super;
  return bfMatIsTransposed(&super->super) ?
    super->rowOffset[bfMatBlockGetNumRowBlocks(super)] :
    super->colOffset[bfMatBlockGetNumColBlocks(super)];
}

BfMatBlockCoo *bfMatBlockCooGetRowRange(BfMatBlockCoo *,
                                        BfSize, BfSize) {
  assert(false);
  return NULL;
}

BfMatBlockCoo *bfMatBlockCooGetColRange(BfMatBlockCoo *,
                                        BfSize, BfSize) {
  assert(false);
  return NULL;
}

void bfMatBlockCooSetRowRange(BfMatBlockCoo *, BfSize, BfSize, BfMat const *) {
  assert(false);
}

void bfMatBlockCooAddInplace(BfMatBlockCoo *, BfMat const *) {
  assert(false);
}

BfMat *bfMatBlockCooMul(BfMatBlockCoo const *op1, BfMat const *op2) {
  BEGIN_ERROR_HANDLING();

  BfSize numRows = bfMatGetNumRows(bfMatBlockCooGetMatConstPtr(op1));
  BfSize numCols = bfMatGetNumCols(op2);
  BfSize numBlocks = bfMatBlockCooNumBlocks(op1);

  BfMat *result = NULL;
  BfMat *block = NULL;
  BfMat *op2Rows = NULL;
  BfMat *tmp = NULL;
  BfMat *resultRows = NULL;

  result = bfMatZerosLike(op2, numRows, numCols);
  HANDLE_ERROR();

  for (BfSize k = 0, i0, i1, j0, j1; k < numBlocks; ++k) {
    i0 = op1->super.rowOffset[op1->rowInd[k]];
    i1 = op1->super.rowOffset[op1->rowInd[k] + 1];
    j0 = op1->super.colOffset[op1->colInd[k]];
    j1 = op1->super.colOffset[op1->colInd[k] + 1];

    block = op1->super.block[k];
    op2Rows = bfMatGetRowRange((BfMat *)op2, j0, j1);
    tmp = bfMatMul(block, op2Rows);
    resultRows = bfMatGetRowRange(result, i0, i1);
    bfMatAddInplace(resultRows, tmp);

    bfMatDeinitAndDelete(&tmp);
  }

  END_ERROR_HANDLING()
    bfMatDeinitAndDelete(&result);

  return result;

}

BfMat *bfMatBlockCooLstSq(BfMatBlockCoo const *, BfMat const *) {
  /* TODO: not sure whether this should be implemented or not... see:
   * https://en.wikipedia.org/wiki/Block_matrix_pseudoinverse */
  assert(false);
  return NULL;
}

BfSize bfMatBlockCooNumBlocks(BfMatBlockCoo const *mat) {
  return mat->numBlocks;
}
