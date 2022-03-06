#include "mat_block_diag.h"

#include <assert.h>
#include <stdlib.h>

#include "error.h"
#include "error_macros.h"

static BfMatVtable matVtbl = {
  .deinit = (__typeof__(&bfMatDeinit))bfMatBlockDiagDeinit,
  .delete = (__typeof__(&bfMatDelete))bfMatBlockDiagDelete,
  .deinitAndDelete = (__typeof__(&bfMatDeinitAndDelete))bfMatBlockDiagDeinitAndDelete,
  .getType = (__typeof__(&bfMatGetType))bfMatBlockDiagGetType,
  .numBytes = (__typeof__(&bfMatNumBytes))bfMatBlockDiagNumBytes,
  .save = (__typeof__(&bfMatSave))bfMatBlockDiagSave,
  .mul = (__typeof__(&bfMatMul))bfMatBlockDiagMul,
  .lstSq = (__typeof__(&bfMatLstSq))bfMatBlockDiagLstSq
};

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

void bfMatBlockDiagDeinit(BfMatBlockDiag *mat) {
  bfMatBlockDeinit(&mat->super);
}

void bfMatBlockDiagDelete(BfMatBlockDiag **mat) {
  free(*mat);
  *mat = NULL;
}

void bfMatBlockDiagDeinitAndDelete(BfMatBlockDiag **mat) {
  bfMatBlockDiagDeinit(*mat);
  bfMatBlockDiagDelete(mat);
}

BfMatType bfMatBlockDiagGetType(BfMatBlockDiag const *mat) {
  return BF_MAT_TYPE_BLOCK_DIAG;
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

BfMat *bfMatBlockDiagMul(BfMatBlockDiag const *op1, BfMat const *op2) {
  (void)op1;
  (void)op2;
  assert(false);
  return NULL;
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
