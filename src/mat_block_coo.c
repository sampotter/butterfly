#include <bf/mat_block_coo.h>

#include <assert.h>
#include <stdlib.h>

#include <bf/error.h>
#include <bf/error_macros.h>

#define INTERFACE BF_INTERFACE_Mat
BF_DEFINE_VTABLE(Mat, MatBlockCoo)
#undef INTERFACE

#define INTERFACE BF_INTERFACE_MatBlock
BF_DEFINE_VTABLE(MatBlock, MatBlockCoo)
#undef INTERFACE

BfMat *bfMatBlockCooToMat(BfMatBlockCoo *matBlockCoo) {
  return &matBlockCoo->super.super;
}

BfMat const *bfMatBlockCooConstToMatConst(BfMatBlockCoo const *matBlockCoo) {
  return &matBlockCoo->super.super;
}

BfMatBlockCoo const *bfMatConstToMatBlockCooConst(BfMat const *mat) {
  if (!bfMatInstanceOf(mat, BF_MAT_TYPE_BLOCK_COO)) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return NULL;
  } else {
    return (BfMatBlockCoo const *)mat;
  }
}

BfMatBlockCoo const *bfMatBlockConstToMatBlockCooConst(BfMatBlock const *matBlock) {
  return bfMatConstToMatBlockCooConst(&matBlock->super);
}

BfMatBlockCoo *bfMatBlockCooNew() {
  return malloc(sizeof(BfMatBlockCoo));
}

void bfMatBlockCooInit(BfMatBlockCoo *mat, BfSize numBlockRows,
                       BfSize numBlockCols, BfSize numBlocks)
{
  BEGIN_ERROR_HANDLING();

  bfMatBlockInit(&mat->super,
                 &MatVtbl, &MatBlockVtbl,
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

void bfMatBlockCooDelete(BfMat **mat) {
  bfMatBlockCooDeinitAndDealloc((BfMatBlockCoo **)mat);
}

BfMat *bfMatBlockCooEmptyLike(BfMat const *mat, BfSize numRows, BfSize numCols) {
  (void)mat; (void)numRows; (void)numCols;
  assert(false);
  return NULL;
}

BfMat *bfMatBlockCooZerosLike(BfMat const *mat, BfSize numRows, BfSize numCols) {
  (void)mat; (void)numRows; (void)numCols;
  assert(false);
  return NULL;
}

void bfMatBlockCooDeinit(BfMatBlockCoo *mat) {
  (void)mat;
  assert(false);
}

void bfMatBlockCooDealloc(BfMatBlockCoo **mat) {
  free(*mat);
  *mat = NULL;
}

void bfMatBlockCooDeinitAndDealloc(BfMatBlockCoo **mat) {
  bfMatBlockCooDeinit(*mat);
  bfMatBlockCooDealloc(mat);
}

BfMatType bfMatBlockCooGetType(BfMat const *mat) {
  (void)mat;
  return BF_MAT_TYPE_BLOCK_COO;
}

bool bfMatBlockCooInstanceOf(BfMat const *mat, BfMatType matType) {
  BfMat const *parent = bfMatBlockCooConstToMatConst((BfMatBlockCoo const *)mat);
  return bfMatTypeDerivedFrom(bfMatGetType(parent), matType);
}

BfSize bfMatBlockCooNumBytes(BfMat const *mat) {
  BfMatBlockCoo const *matBlockCoo = bfMatConstToMatBlockCooConst(mat);

  BfSize num_bytes = 0;

  /* memory occupied by mat itself */
  num_bytes += sizeof(BfMatBlockCoo);

  /* memory occupied by blocks */
  for (BfSize i = 0; i < matBlockCoo->numBlocks; ++i)
    num_bytes += bfMatNumBytes(matBlockCoo->super.block[i]);

  /* mat->super.rowOffset and mat->super.colOffset */
  num_bytes += (mat->numRows + 1)*sizeof(BfSize);
  num_bytes += (mat->numCols + 1)*sizeof(BfSize);

  /* mat->rowInd and mat->colInd */
  num_bytes += 2*matBlockCoo->numBlocks*sizeof(BfSize);

  return num_bytes;
}

void bfMatBlockCooSave(BfMat const *mat, char const *path) {
  (void)mat; (void)path;
  assert(false);
}

BfSize bfMatBlockCooGetNumRows(BfMat const *mat) {
  BfMatBlock const *matBlock = bfMatConstToMatBlockConst(mat);
  return bfMatIsTransposed(mat) ?
    matBlock->colOffset[bfMatBlockGetNumColBlocks(matBlock)] :
    matBlock->rowOffset[bfMatBlockGetNumRowBlocks(matBlock)];
}

BfSize bfMatBlockCooGetNumCols(BfMat const *mat) {
  BfMatBlock const *matBlock = bfMatConstToMatBlockConst(mat);
  return bfMatIsTransposed(mat) ?
    matBlock->rowOffset[bfMatBlockGetNumRowBlocks(matBlock)] :
    matBlock->colOffset[bfMatBlockGetNumColBlocks(matBlock)];
}

BfMat *bfMatBlockCooGetRowRange(BfMat *mat, BfSize i0, BfSize i1) {
  (void)mat; (void)i0; (void)i1;
  assert(false);
  return NULL;
}

BfMat *bfMatBlockCooGetColRange(BfMat *mat, BfSize j0, BfSize j1) {
  (void)mat; (void)j0; (void)j1;
  assert(false);
  return NULL;
}

void bfMatBlockCooSetRowRange(BfMat *mat, BfSize i0, BfSize i1, BfMat const *otherMat) {
  (void)mat; (void)i0; (void)i1; (void)otherMat;
  assert(false);
}

void bfMatBlockCooAddInplace(BfMat *mat, BfMat const *otherMat) {
  (void)mat;
  (void)otherMat;
  assert(false);
}

BfMat *bfMatBlockCooMul(BfMat const *op1, BfMat const *op2) {
  BEGIN_ERROR_HANDLING();

  BfMatBlock const *matBlock1 = bfMatConstToMatBlockConst(op1);
  BfMatBlockCoo const *matBlockCoo1 = bfMatConstToMatBlockCooConst(op1);

  BfSize numRows = bfMatGetNumRows(op1);
  BfSize numCols = bfMatGetNumCols(op2);
  BfSize numBlocks = bfMatBlockCooNumBlocks(matBlock1);

  BfMat *result = NULL;
  BfMat *block = NULL;
  BfMat *op2Rows = NULL;
  BfMat *tmp = NULL;
  BfMat *resultRows = NULL;

  result = bfMatZerosLike(op2, numRows, numCols);
  HANDLE_ERROR();

  for (BfSize k = 0, i0, i1, j0, j1; k < numBlocks; ++k) {
    i0 = matBlock1->rowOffset[matBlockCoo1->rowInd[k]];
    i1 = matBlock1->rowOffset[matBlockCoo1->rowInd[k] + 1];
    j0 = matBlock1->colOffset[matBlockCoo1->colInd[k]];
    j1 = matBlock1->colOffset[matBlockCoo1->colInd[k] + 1];

    block = matBlock1->block[k];
    op2Rows = bfMatGetRowRange((BfMat *)op2, j0, j1);
    tmp = bfMatMul(block, op2Rows);
    resultRows = bfMatGetRowRange(result, i0, i1);
    bfMatAddInplace(resultRows, tmp);

    bfMatDelete(&tmp);
  }

  END_ERROR_HANDLING()
    bfMatDelete(&result);

  return result;

}

BfMat *bfMatBlockCooLstSq(BfMat const *mat, BfMat const *otherMat) {
  /* TODO: not sure whether this should be implemented or not... see:
   * https://en.wikipedia.org/wiki/Block_matrix_pseudoinverse */
  (void)mat; (void)otherMat;
  assert(false);
  return NULL;
}

BfSize bfMatBlockCooNumBlocks(BfMatBlock const *mat) {
  return bfMatBlockConstToMatBlockCooConst(mat)->numBlocks;
}

BfSize bfMatBlockCooGetNumRowBlocks(BfMatBlock const *matBlock) {
  (void)matBlock;
  assert(false);
}

BfSize bfMatBlockCooGetNumColBlocks(BfMatBlock const *matBlock) {
  (void)matBlock;
  assert(false);
}

BfSize bfMatBlockCooGetNumBlockRows(BfMatBlock const *matBlock, BfSize i) {
  (void)matBlock;
  (void)i;
  assert(false);
}

BfSize bfMatBlockCooGetNumBlockCols(BfMatBlock const *matBlock, BfSize j) {
  (void)matBlock;
  (void)j;
  assert(false);
}
