#include <bf/mat_block.h>

#include <stdlib.h>
#include <string.h>

#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/util.h>

/** Interface: MatBlock */

BfMat *bfMatBlockGetBlockCopy(BfMatBlock const *matBlock, BfSize i, BfSize j) {
  return matBlock->vtbl->GetBlockCopy(matBlock, i, j);
}

/** Upcasting: MatBlock -> Mat */

BfMat const *bfMatBlockConstToMatConst(BfMatBlock const *matBlock) {
  return &matBlock->super;
}

/** Downcasting: Mat -> MatBlock */

BfMatBlock *bfMatToMatBlock(BfMat *mat) {
  if (!bfMatInstanceOf(mat, BF_TYPE_MAT_BLOCK)) {
    bfSetError(BF_ERROR_RUNTIME_ERROR);
    return NULL;
  } else {
    return (BfMatBlock *)mat;
  }
}

BfMatBlock const *bfMatConstToMatBlockConst(BfMat const *mat) {
  if (!bfMatInstanceOf(mat, BF_TYPE_MAT_BLOCK)) {
    bfSetError(BF_ERROR_RUNTIME_ERROR);
    return NULL;
  } else {
    return (BfMatBlock const *)mat;
  }
}

/** Implementation: MatBlock */

void bfMatBlockInvalidate(BfMatBlock *matBlock) {
  bfMatInvalidate(&matBlock->super);

  matBlock->vtbl = NULL;
  matBlock->block = NULL;
  matBlock->rowOffset = NULL;
  matBlock->colOffset = NULL;
}

void bfMatBlockInit(BfMatBlock *mat,
                    BfMatVtable *matVtbl, BfMatBlockVtable *matBlockVtbl,
                    BfSize numBlocks, BfSize numBlockRows, BfSize numBlockCols)
{
  BEGIN_ERROR_HANDLING();

  bfMatInit(&mat->super, matVtbl, numBlockRows, numBlockCols);
  HANDLE_ERROR();

  mat->vtbl = matBlockVtbl;

  mat->block = malloc(numBlocks*sizeof(BfMat *));
  if (mat->block == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);
  memset(mat->block, 0x0, numBlocks*sizeof(BfMat *));

  mat->rowOffset = malloc((numBlockRows + 1)*sizeof(BfSize));
  if (mat->rowOffset == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);
  bfSizeSetConstant(numBlockRows + 1, mat->rowOffset, BF_SIZE_BAD_VALUE);

  mat->colOffset = malloc((numBlockCols + 1)*sizeof(BfSize));
  if (mat->colOffset == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);
  bfSizeSetConstant(numBlockCols + 1, mat->colOffset, BF_SIZE_BAD_VALUE);

  END_ERROR_HANDLING()
    bfMatBlockDeinit(mat);
}

void bfMatBlockDelete(BfMat **mat) {
  bfMatBlockDeinitAndDealloc((BfMatBlock **)mat);
}

BfSize bfMatBlockGetNumBlockRows(BfMatBlock const *mat, BfSize i) {
  return mat->rowOffset[i + 1] - mat->rowOffset[i];
}

BfSize bfMatBlockGetNumBlockCols(BfMatBlock const *mat, BfSize j) {
  return mat->colOffset[j + 1] - mat->colOffset[j];
}

BfSize bfMatBlockNumBlocks(BfMatBlock const *mat) {
  return mat->vtbl->NumBlocks(mat);
}

BfSize bfMatBlockGetNumRowBlocks(BfMatBlock const *mat) {
  return mat->super.numRows;
}

BfSize bfMatBlockGetNumColBlocks(BfMatBlock const *mat) {
  return mat->super.numCols;
}

BfSize bfMatBlockFindRowBlock(BfMatBlock const *mat, BfSize i0) {
  BfSize p = 0;
  while (p < mat->super.numRows && mat->rowOffset[p] < i0)
    ++p;
  return mat->rowOffset[p] <= i0 && i0 < mat->rowOffset[p + 1] ?
    p :
    BF_SIZE_BAD_VALUE;
}

void bfMatBlockDeinit(BfMatBlock *mat) {
  bfMatDeinit(&mat->super);

  free(mat->block);
  free(mat->rowOffset);
  free(mat->colOffset);

#if BF_DEBUG
  mat->vtbl = NULL;
  mat->block = NULL;
  mat->rowOffset = NULL;
  mat->colOffset = NULL;
#endif
}

void bfMatBlockDealloc(BfMatBlock **mat) {
  free(*mat);
  *mat = NULL;
}

void bfMatBlockDeinitAndDealloc(BfMatBlock **mat) {
  bfMatBlockDeinit(*mat);
  bfMatBlockDealloc(mat);
}

void bfMatBlockDelete(BfMat **mat) {
  bfMatBlockDeinitAndDealloc((BfMatBlock **)mat);
}
