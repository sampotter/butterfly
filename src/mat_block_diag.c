#include <bf/mat_block_diag.h>

#include <assert.h>
#include <stdlib.h>

#include <bf/error.h>
#include <bf/error_macros.h>

BF_DEFINE_MAT_VTABLE(MatBlockDiag);

static BfMatBlockVtable matBlockVtbl = {
  .numBlocks = (__typeof__(&bfMatBlockNumBlocks))bfMatBlockDiagNumBlocks
};

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

  bfMatBlockInit(&mat->super, &matVtbl, &matBlockVtbl, numBlocks, numRows, numCols);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {}
}

BfMat *bfMatBlockDiagGetMatPtr(BfMatBlockDiag *mat) {
  return &mat->super.super;
}

BfMat const *bfMatBlockDiagGetMatConstPtr(BfMatBlockDiag const *mat) {
  return &mat->super.super;
}

void bfMatBlockDiagDelete(BfMatBlockDiag **mat) {
  bfMatBlockDiagDeinitAndDealloc(mat);
}

BfMatBlockDiag *bfMatBlockDiagEmptyLike(BfMatBlockDiag const *, BfSize, BfSize) {
  assert(false);
  return NULL;
}

BfMatBlockDiag *bfMatBlockDiagZerosLike(BfMatBlockDiag const *, BfSize, BfSize) {
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

BfMatType bfMatBlockDiagGetType(BfMatBlockDiag const *mat) {
  return BF_MAT_TYPE_BLOCK_DIAG;
}

bool bfMatBlockDiagInstanceOf(BfMatBlockDiag const *mat, BfMatType matType) {
  BfMat const *parent = bfMatBlockDiagGetMatConstPtr(mat);
  return bfMatTypeDerivedFrom(bfMatGetType(parent), matType);
}

BfSize bfMatBlockDiagNumBytes(BfMatBlockDiag const *mat) {
  (void)mat;
  assert(false);
  return 0;
}

void bfMatBlockDiagSave(BfMatBlockDiag const *mat, char const *path) {
  (void)mat;
  (void)path;
  assert(false);
}

BfSize bfMatBlockDiagGetNumRows(BfMatBlockDiag const *mat) {
  BfMatBlock const *super = &mat->super;
  return bfMatIsTransposed(&super->super) ?
    super->colOffset[bfMatBlockGetNumColBlocks(super)] :
    super->rowOffset[bfMatBlockGetNumRowBlocks(super)];
}

BfSize bfMatBlockDiagGetNumCols(BfMatBlockDiag const *mat) {
  BfMatBlock const *super = &mat->super;
  return bfMatIsTransposed(&super->super) ?
    super->rowOffset[bfMatBlockGetNumRowBlocks(super)] :
    super->colOffset[bfMatBlockGetNumColBlocks(super)];
}

BfMatBlockDiag *bfMatBlockDiagGetRowRange(BfMatBlockDiag *, BfSize, BfSize) {
  assert(false);
  return NULL;
}

BfMatBlockDiag *bfMatBlockDiagGetColRange(BfMatBlockDiag *, BfSize, BfSize) {
  assert(false);
  return NULL;
}

void bfMatBlockDiagSetRowRange(BfMatBlockDiag *, BfSize, BfSize, BfMat const *) {
  assert(false);
}

void bfMatBlockDiagAddInplace(BfMatBlockDiag *, BfMat const *) {
  assert(false);
}

BfMat *bfMatBlockDiagMul(BfMatBlockDiag const *op1, BfMat const *op2) {
  BEGIN_ERROR_HANDLING();

  BfSize numRows = bfMatGetNumRows(bfMatBlockDiagGetMatConstPtr(op1));
  BfSize numCols = bfMatGetNumCols(op2);
  BfSize numBlocks = bfMatBlockDiagNumBlocks(op1);

  BfMat *result = NULL;
  BfMat *block = NULL;
  BfMat *op2Rows = NULL;
  BfMat *resultRows = NULL;

  result = bfMatEmptyLike(op2, numRows, numCols);
  HANDLE_ERROR();

  for (BfSize i = 0, i0, i1, j0, j1; i < numBlocks; ++i) {
    i0 = op1->super.rowOffset[i];
    i1 = op1->super.rowOffset[i + 1];
    j0 = op1->super.colOffset[i];
    j1 = op1->super.colOffset[i + 1];

    block = op1->super.block[i];
    op2Rows = bfMatGetRowRange((BfMat *)op2, j0, j1);
    resultRows = bfMatMul(block, op2Rows);
    bfMatSetRowRange(result, i0, i1, resultRows);

    bfMatDelete(&resultRows);
  }

  END_ERROR_HANDLING()
    bfMatDelete(&result);

  return result;
}

BfMat *bfMatBlockDiagLstSq(BfMatBlockDiag const *op1, BfMat const *op2) {
  (void)op1;
  (void)op2;
  assert(false);
  return NULL;
}

BfSize bfMatBlockDiagNumBlocks(BfMatBlockDiag const *mat) {
  BfSize m = mat->super.super.numRows;
  BfSize n = mat->super.super.numCols;
  return m < n ? m : n;
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
