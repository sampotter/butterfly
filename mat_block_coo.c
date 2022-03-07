#include "mat_block_coo.h"

#include <assert.h>
#include <stdlib.h>

#include "error.h"
#include "error_macros.h"

static BfMatVtable matVtbl = {
  .deinit = (__typeof__(&bfMatDeinit))bfMatBlockCooDeinit,
  .delete = (__typeof__(&bfMatDelete))bfMatBlockCooDelete,
  .deinitAndDelete = (__typeof__(&bfMatDeinitAndDelete))bfMatBlockCooDeinitAndDelete,
  .getType = (__typeof__(&bfMatGetType))bfMatBlockCooGetType,
  .numBytes = (__typeof__(&bfMatNumBytes))bfMatBlockCooNumBytes,
  .save = (__typeof__(&bfMatSave))bfMatBlockCooSave,
  .mul = (__typeof__(&bfMatMul))bfMatBlockCooMul,
  .lstSq = (__typeof__(&bfMatLstSq))bfMatBlockCooLstSq
};

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

BfMatType bfMatBlockCooGetType(BfMatBlockCoo *mat) {
  (void)mat;
  return BF_MAT_TYPE_BLOCK_COO;
}

BfSize bfMatBlockCooNumBytes(BfMatBlockCoo *mat) {
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

void bfMatBlockCooSave(BfMatBlockCoo const *mat, char const *path) {
  (void)mat;
  (void)path;
  assert(false);
}

BfMat *bfMatBlockCooMul(BfMatBlockCoo const *op1, BfMat const *op2) {
  (void)op1;
  (void)op2;
  assert(false);
  return NULL;
}

BfMat *bfMatBlockCooLstSq(BfMatBlockCoo const *lhs, BfMat const *rhs) {
  (void)lhs;
  (void)rhs;
  assert(false);
  return NULL;

  // TODO: not sure whether this should be implemented or not... see:
  //   https://en.wikipedia.org/wiki/Block_matrix_pseudoinverse
}

BfSize bfMatBlockCooNumBlocks(BfMatBlockCoo const *mat) {
  return mat->numBlocks;
}
