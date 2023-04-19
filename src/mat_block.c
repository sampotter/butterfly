#include <bf/mat_block.h>

#include <bf/assert.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mem.h>
#include <bf/util.h>

/** Interface: MatBlock */

BfSize bfMatBlockNumBlocks(BfMatBlock const *mat) {
  return mat->vtbl->NumBlocks(mat);
}

BfSize bfMatBlockGetNumRowBlocks(BfMatBlock const *matBlock) {
  return matBlock->vtbl->GetNumRowBlocks(matBlock);
}

BfSize bfMatBlockGetNumColBlocks(BfMatBlock const *matBlock) {
  return matBlock->vtbl->GetNumColBlocks(matBlock);
}

BfSize bfMatBlockGetNumBlockRows(BfMatBlock const *mat, BfSize i) {
  BF_ASSERT(mat->vtbl->GetNumBlockRows == NULL); // TODO: safeguard
  return mat->rowOffset[i + 1] - mat->rowOffset[i];
}

BfSize bfMatBlockGetNumBlockCols(BfMatBlock const *mat, BfSize j) {
  BF_ASSERT(mat->vtbl->GetNumBlockCols == NULL); // TODO: safeguard
  return mat->colOffset[j + 1] - mat->colOffset[j];
}

BfSize bfMatBlockGetRowOffset(BfMatBlock const *matBlock, BfSize i) {
  return matBlock->vtbl->GetRowOffset(matBlock, i);
}

BfSize bfMatBlockGetColOffset(BfMatBlock const *matBlock, BfSize j) {
  return matBlock->vtbl->GetColOffset(matBlock, j);
}

BfMat const *bfMatBlockGetBlockConst(BfMatBlock const *matBlock, BfSize i, BfSize j) {
  return matBlock->vtbl->GetBlockConst(matBlock, i, j);
}

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

  mat->block = bfMemAllocAndZero(numBlocks, sizeof(BfMat *));
  if (mat->block == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  mat->rowOffset = bfMemAlloc(numBlockRows + 1, sizeof(BfSize));
  if (mat->rowOffset == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);
  bfSizeSetConstant(numBlockRows + 1, mat->rowOffset, BF_SIZE_BAD_VALUE);

  mat->colOffset = bfMemAlloc(numBlockCols + 1, sizeof(BfSize));
  if (mat->colOffset == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);
  bfSizeSetConstant(numBlockCols + 1, mat->colOffset, BF_SIZE_BAD_VALUE);

  END_ERROR_HANDLING()
    bfMatBlockDeinit(mat);
}

void bfMatBlockDeinit(BfMatBlock *mat) {
  bfMatDeinit(&mat->super);

  bfMemFree(mat->block);
  bfMemFree(mat->rowOffset);
  bfMemFree(mat->colOffset);

#if BF_DEBUG
  mat->vtbl = NULL;
  mat->block = NULL;
  mat->rowOffset = NULL;
  mat->colOffset = NULL;
#endif
}

void bfMatBlockDealloc(BfMatBlock **mat) {
  bfMemFree(*mat);
  *mat = NULL;
}

void bfMatBlockDeinitAndDealloc(BfMatBlock **mat) {
  bfMatBlockDeinit(*mat);
  bfMatBlockDealloc(mat);
}

void bfMatBlockDelete(BfMat **mat) {
  bfMatBlockDeinitAndDealloc((BfMatBlock **)mat);
}

BfSize bfMatBlockFindRowBlock(BfMatBlock const *mat, BfSize i0) {
  BfSize p = 0;
  while (p < mat->super.numRows && mat->rowOffset[p] < i0)
    ++p;
  return mat->rowOffset[p] <= i0 && i0 < mat->rowOffset[p + 1] ?
    p :
    BF_SIZE_BAD_VALUE;
}

BfSize getMaxDepthRec(BfMatBlock const *matBlock, BfSize maxDepth) {
  BfSize newMaxDepth = maxDepth;

  BfSize numBlocks = bfMatBlockNumBlocks(matBlock);
  for (BfSize i = 0; i < numBlocks; ++i) {
    BfSize depth = maxDepth + 1;
    if (bfMatIsBlock(matBlock->block[i])) {
      BfMatBlock const *childBlock = bfMatConstToMatBlockConst(matBlock->block[i]);
      depth = getMaxDepthRec(childBlock, maxDepth + 1);
    }
    if (depth > newMaxDepth) newMaxDepth = depth;
  }

  return newMaxDepth;
}

BfSize bfMatBlockGetMaxDepth(BfMatBlock const *matBlock) {
  return getMaxDepthRec(matBlock, 0);
}
