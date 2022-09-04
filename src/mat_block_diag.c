#include <bf/mat_block_diag.h>

#include <assert.h>
#include <stdlib.h>

#include <bf/error.h>
#include <bf/error_macros.h>

#define INTERFACE BF_INTERFACE_Mat
BF_DEFINE_VTABLE(Mat, MatBlockDiag)
#undef INTERFACE

#define INTERFACE BF_INTERFACE_MatBlock
BF_DEFINE_VTABLE(MatBlock, MatBlockDiag)
#undef INTERFACE

BfMat *bfMatBlockDiagToMat(BfMatBlockDiag *matBlockDiag) {
  return &matBlockDiag->super.super;
}

BfMatBlockDiag const *bfMatConstToMatBlockDiagConst(BfMat const *mat) {
  if (!bfMatInstanceOf(mat, BF_MAT_TYPE_BLOCK_DIAG)) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return NULL;
  } else {
    return (BfMatBlockDiag const *)mat;
  }
}

BfMatBlockDiag *bfMatBlockDiagNew() {
  BEGIN_ERROR_HANDLING();

  BfMatBlockDiag *mat = malloc(sizeof(BfMatBlockDiag));
  if (mat == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING() {}

  return mat;
}

void bfMatBlockDiagInit(BfMatBlockDiag *mat, BfSize numRows, BfSize numCols) {
  BEGIN_ERROR_HANDLING();

  BfSize numBlocks = numRows < numCols ? numRows : numCols;

  bfMatBlockInit(&mat->super, &MatVtbl, &MatBlockVtbl, numBlocks, numRows, numCols);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {}
}

void bfMatBlockDiagDelete(BfMat **mat) {
  bfMatBlockDiagDeinitAndDealloc((BfMatBlockDiag **)mat);
}

BfMat *bfMatBlockDiagEmptyLike(BfMat const *mat, BfSize numRows, BfSize numCols) {
  (void)mat;
  (void)numRows;
  (void)numCols;
  assert(false);
  return NULL;
}

BfMat *bfMatBlockDiagZerosLike(BfMat const *mat, BfSize numRows, BfSize numCols) {
  (void)mat;
  (void)numRows;
  (void)numCols;
  assert(false);
  return NULL;
}

void bfMatBlockDiagDeinit(BfMatBlockDiag *mat) {
  bfMatBlockDeinit(&mat->super);
}

void bfMatBlockDiagDealloc(BfMatBlockDiag **mat) {
  free(*mat);
  *mat = NULL;
}

void bfMatBlockDiagDeinitAndDealloc(BfMatBlockDiag **mat) {
  bfMatBlockDiagDeinit(*mat);
  bfMatBlockDiagDealloc(mat);
}

BfMatType bfMatBlockDiagGetType(BfMat const *mat) {
  (void)mat;
  return BF_MAT_TYPE_BLOCK_DIAG;
}

bool bfMatBlockDiagInstanceOf(BfMat const *mat, BfMatType matType) {
  return bfMatTypeDerivedFrom(bfMatGetType(mat), matType);
}

BfSize bfMatBlockDiagNumBytes(BfMat const *mat) {
  (void)mat;
  assert(false);
  return 0;
}

void bfMatBlockDiagSave(BfMat const *mat, char const *path) {
  (void)mat;
  (void)path;
  assert(false);
}

BfSize bfMatBlockDiagGetNumRows(BfMat const *mat) {
  BfMatBlock const *matBlock = bfMatConstToMatBlockConst(mat);
  return bfMatIsTransposed(mat) ?
    matBlock->colOffset[bfMatBlockGetNumColBlocks(matBlock)] :
    matBlock->rowOffset[bfMatBlockGetNumRowBlocks(matBlock)];
}

BfSize bfMatBlockDiagGetNumCols(BfMat const *mat) {
  BfMatBlock const *matBlock = bfMatConstToMatBlockConst(mat);
  return bfMatIsTransposed(mat) ?
    matBlock->rowOffset[bfMatBlockGetNumRowBlocks(matBlock)] :
    matBlock->colOffset[bfMatBlockGetNumColBlocks(matBlock)];
}

BfMat *bfMatBlockDiagGetRowRange(BfMat *mat, BfSize i0, BfSize i1) {
  (void)mat; (void)i0; (void)i1;
  assert(false);
  return NULL;
}

BfMat *bfMatBlockDiagGetColRange(BfMat *mat, BfSize j0, BfSize j1) {
  (void)mat; (void)j0; (void)j1;
  assert(false);
  return NULL;
}

void bfMatBlockDiagSetRowRange(BfMat *mat, BfSize i0, BfSize i1, BfMat const *otherMat) {
  (void)mat; (void)i0; (void)i1; (void)otherMat;
  assert(false);
}

void bfMatBlockDiagAddInplace(BfMat *mat, BfMat const *otherMat) {
  (void)mat; (void)otherMat;
  assert(false);
}

BfMat *bfMatBlockDiagMul(BfMat const *mat, BfMat const *other) {
  BEGIN_ERROR_HANDLING();

  BfMatBlock const *matBlock = bfMatConstToMatBlockConst(mat);

  BfSize numRows = bfMatGetNumRows(mat);
  BfSize numCols = bfMatGetNumCols(other);
  BfSize numBlocks = bfMatBlockDiagNumBlocks(matBlock);

  BfMat *result = NULL;
  BfMat *block = NULL;
  BfMat *op2Rows = NULL;
  BfMat *resultRows = NULL;

  result = bfMatEmptyLike(other, numRows, numCols);
  HANDLE_ERROR();

  for (BfSize i = 0, i0, i1, j0, j1; i < numBlocks; ++i) {
    i0 = matBlock->rowOffset[i];
    i1 = matBlock->rowOffset[i + 1];
    j0 = matBlock->colOffset[i];
    j1 = matBlock->colOffset[i + 1];

    block = matBlock->block[i];
    op2Rows = bfMatGetRowRange((BfMat *)other, j0, j1);
    resultRows = bfMatMul(block, op2Rows);
    bfMatSetRowRange(result, i0, i1, resultRows);

    bfMatDelete(&resultRows);
  }

  END_ERROR_HANDLING()
    bfMatDelete(&result);

  return result;
}

BfMat *bfMatBlockDiagLstSq(BfMat const *op1, BfMat const *op2) {
  (void)op1;
  (void)op2;
  assert(false);
  return NULL;
}

BfSize bfMatBlockDiagNumBlocks(BfMatBlock const *matBlock) {
  BfMat const *mat = bfMatBlockConstToMatConst(matBlock);
  return mat->numRows < mat->numCols ? mat->numRows : mat->numCols;
}

BfMat const *
bfMatBlockDiagGetBlock(BfMatBlockDiag const *mat, BfSize i, BfSize j) {
  BEGIN_ERROR_HANDLING();

  BfMat const *block = NULL;

  if (i >= bfMatBlockGetNumRowBlocks(&mat->super))
    RAISE_ERROR(BF_ERROR_OUT_OF_RANGE);

  if (j >= bfMatBlockGetNumColBlocks(&mat->super))
    RAISE_ERROR(BF_ERROR_OUT_OF_RANGE);

  if (i != j)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  block = mat->super.block[i];

  END_ERROR_HANDLING()
    block = NULL;

  return block;
}

BfSize bfMatBlockDiagGetNumRowBlocks(BfMatBlock const *mat) {
  (void)mat;
  assert(false);
}

BfSize bfMatBlockDiagGetNumColBlocks(BfMatBlock const *mat) {
  (void)mat;
  assert(false);
}

BfSize bfMatBlockDiagGetNumBlockRows(BfMatBlock const *mat, BfSize i) {
  (void)mat;
  (void)i;
  assert(false);
}

BfSize bfMatBlockDiagGetNumBlockCols(BfMatBlock const *mat, BfSize j) {
  (void)mat;
  (void)j;
  assert(false);
}
