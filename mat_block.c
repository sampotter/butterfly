#include "mat_block.h"

#include <stdlib.h>
#include <string.h>

#include "error.h"
#include "error_macros.h"
#include "util.h"

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

BfSize bfMatBlockGetNumBlockRows(BfMatBlock const *mat, BfSize i) {
  return mat->rowOffset[i + 1] - mat->rowOffset[i];
}

BfSize bfMatBlockGetNumBlockCols(BfMatBlock const *mat, BfSize j) {
  return mat->colOffset[j + 1] - mat->colOffset[j];
}

BfSize bfMatBlockNumBlocks(BfMatBlock const *mat) {
  return mat->vtbl->numBlocks(mat);
}

void bfMatBlockDeinit(BfMatBlock *mat) {
  bfMatDeinit(&mat->super);

  mat->vtbl = NULL;

  free(mat->block);
  mat->block = NULL;

  free(mat->rowOffset);
  mat->rowOffset = NULL;

  free(mat->colOffset);
  mat->colOffset = NULL;
}

void bfMatBlockDelete(BfMatBlock **mat) {
  free(*mat);
  *mat = NULL;
}

void bfMatBlockDeinitAndDelete(BfMatBlock **mat) {
  bfMatBlockDeinit(*mat);
  bfMatBlockDelete(mat);
}
