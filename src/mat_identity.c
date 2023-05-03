#include <bf/mat_identity.h>

#include <bf/assert.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mem.h>

/** Interface: Mat */

static BfMatVtable MAT_VTABLE = {
  .Copy = (__typeof__(&bfMatCopy))bfMatIdentityCopy,
  .Steal = (__typeof__(&bfMatSteal))bfMatIdentitySteal,
  .Delete = (__typeof__(&bfMatDelete))bfMatIdentityDelete,
  .GetType = (__typeof__(&bfMatGetType))bfMatIdentityGetType,
  .NumBytes = (__typeof__(&bfMatNumBytes))bfMatIdentityNumBytes,
  .GetNumRows = (__typeof__(&bfMatGetNumRows))bfMatIdentityGetNumRows,
  .GetNumCols = (__typeof__(&bfMatGetNumCols))bfMatIdentityGetNumCols,
  .GetRowRangeCopy = (__typeof__(&bfMatGetRowRangeCopy))bfMatIdentityGetRowRangeCopy,
  .MulVec = (__typeof__(&bfMatMulVec))bfMatIdentityMulVec,
  .PrintBlocksDeep = (__typeof__(&bfMatPrintBlocksDeep))bfMatIdentityPrintBlocksDeep,
};

BfMat *bfMatIdentityCopy(BfMatIdentity const *matIdentity) {
  BEGIN_ERROR_HANDLING();

  BfSize numRows = bfMatIdentityGetNumRows(matIdentity);
  BfSize numCols = bfMatIdentityGetNumCols(matIdentity);

  if (numRows != numCols)
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);

  BfMatIdentity *copy = bfMatIdentityNew();
  HANDLE_ERROR();

  bfMatIdentityInit(copy, numRows);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {
    BF_ASSERT(false);
  }

  return bfMatIdentityToMat(copy);
}

BfMat *bfMatIdentitySteal(BfMatIdentity *matIdentity) {
  BEGIN_ERROR_HANDLING();

  BfMat *mat = bfMatIdentityToMat(matIdentity);

  if (bfMatIsView(mat))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfMatIdentity *matIdentityNew = bfMatIdentityNew();
  HANDLE_ERROR();

  *matIdentityNew = *matIdentity;

  mat->props |= BF_MAT_PROPS_VIEW;

  END_ERROR_HANDLING() {
    BF_DIE();
  }

  return bfMatIdentityToMat(matIdentityNew);
}

void bfMatIdentityDelete(BfMatIdentity **matIdentity) {
  bfMatIdentityDeinitAndDealloc(matIdentity);
}

BfType bfMatIdentityGetType(BfMatIdentity const *matIdentity) {
  (void)matIdentity;
  return BF_TYPE_MAT_IDENTITY;
}

BfSize bfMatIdentityNumBytes(BfMatIdentity const *matIdentity) {
  (void)matIdentity;
  return 0;
}

BfSize bfMatIdentityGetNumRows(BfMatIdentity const *matIdentity) {
  return matIdentity->super.numRows;
}

BfSize bfMatIdentityGetNumCols(BfMatIdentity const *matIdentity) {
  return matIdentity->super.numCols;
}

BfMat *bfMatIdentityGetRowRangeCopy(BfMatIdentity const *matIdentity, BfSize i0, BfSize i1) {
  BEGIN_ERROR_HANDLING();

  BfMat const *mat = bfMatIdentityConstToMatConst(matIdentity);

  BfMat *rowRangeCopy = NULL;

  BfSize m = bfMatGetNumRows(mat);

  if (i0 == 0 && i1 == m) {
    BfMatIdentity *matIdentity_ = bfMatIdentityNew();
    HANDLE_ERROR();

    bfMatIdentityInit(matIdentity_, m);
    HANDLE_ERROR();

    rowRangeCopy = bfMatIdentityToMat(matIdentity_);
  } else {
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);
  }

  END_ERROR_HANDLING() {}

  return rowRangeCopy;
}

BfVec *bfMatIdentityMulVec(BfMatIdentity const *matIdentity, BfVec const *vec) {
  BEGIN_ERROR_HANDLING();

  BfMat const *mat = bfMatIdentityConstToMatConst(matIdentity);

  if (bfMatGetNumRows(mat) != bfMatGetNumCols(mat))
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);

  BfVec *result = bfVecCopy(vec);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {
    BF_ASSERT(false);
  }

  return result;
}

void bfMatIdentityPrintBlocksDeep(BfMatIdentity const *matIdentity, FILE *fp, BfSize i0, BfSize j0, BfSize depth) {
  BfMat const *mat = bfMatIdentityConstToMatConst(matIdentity);

  BfSize i1 = i0 + bfMatGetNumRows(mat);
  BfSize j1 = j0 + bfMatGetNumCols(mat);

  fprintf(fp, "%u %lu %lu %lu %lu %lu\n", BF_TYPE_MAT_IDENTITY, i0, i1, j0, j1, depth);
}

/** Upcasting: MatIdentity -> Mat */

BfMat *bfMatIdentityToMat(BfMatIdentity *matIdentity) {
  return &matIdentity->super;
}

BfMat const *bfMatIdentityConstToMatConst(BfMatIdentity const *matIdentity) {
  return &matIdentity->super;
}

/** Downcasting: Mat -> MatIdentity */

BfMatIdentity *bfMatToMatIdentity(BfMat *mat) {
  if (!bfMatInstanceOf(mat, BF_TYPE_MAT_IDENTITY)) {
    bfSetError(BF_ERROR_RUNTIME_ERROR);
    return NULL;
  } else {
    return (BfMatIdentity *)mat;
  }
}

/** Implementation: MatIdentity */

BfMatIdentity *bfMatIdentityNew() {
  BEGIN_ERROR_HANDLING();

  BfMatIdentity *matIdentity = bfMemAlloc(1, sizeof(BfMatIdentity));
  if (matIdentity == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING() {}

  return matIdentity;
}

void bfMatIdentityInit(BfMatIdentity *matIdentity, BfSize n) {
  BEGIN_ERROR_HANDLING();

  bfMatInit(&matIdentity->super, &MAT_VTABLE, n, n);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {
    bfMatIdentityDeinit(matIdentity);
  }
}

void bfMatIdentityDeinit(BfMatIdentity *matIdentity) {
  (void)matIdentity;
  BF_ASSERT(false);
}

void bfMatIdentityDealloc(BfMatIdentity **matIdentity) {
  bfMemFree(*matIdentity);
  *matIdentity = NULL;
}

void bfMatIdentityDeinitAndDealloc(BfMatIdentity **matIdentity) {
  bfMatIdentityDeinit(*matIdentity);
  bfMatIdentityDealloc(matIdentity);
}
