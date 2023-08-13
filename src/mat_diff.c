#include <bf/mat_diff.h>

#include <bf/assert.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mem.h>

/** Interface: MatDiff */

static BfMatVtable MAT_VTABLE = {
  .GetType = (__typeof__(&bfMatGetType))bfMatDiffGetType,
  .GetNumRows = (__typeof__(&bfMatGetNumRows))bfMatDiffGetNumRows,
  .GetNumCols = (__typeof__(&bfMatGetNumCols))bfMatDiffGetNumCols,
  .Mul = (__typeof__(&bfMatMul))bfMatDiffMul,
  .ToType = (__typeof__(&bfMatToType))bfMatDiffToType,
};

BfType bfMatDiffGetType(BfMatDiff const *matDiff) {
  (void)matDiff;
  return BF_TYPE_MAT_DIFF;
}

BfSize bfMatDiffGetNumRows(BfMatDiff const *matDiff) {
  BfMat const *mat = bfMatDiffConstToMatConst(matDiff);
  return mat->numRows;
}

BfSize bfMatDiffGetNumCols(BfMatDiff const *matDiff) {
  BfMat const *mat = bfMatDiffConstToMatConst(matDiff);
  return mat->numCols;
}

BfMat *bfMatDiffMul(BfMatDiff const *matDiff, BfMat const *otherMat) {
  BF_ERROR_BEGIN();

  if (bfMatDiffGetNumCols(matDiff) != bfMatGetNumRows(otherMat))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfMat *firstProduct = bfMatMul(matDiff->first, otherMat);
  HANDLE_ERROR();

  BfMat *secondProduct = bfMatMul(matDiff->second, otherMat);
  HANDLE_ERROR();

  BfMat *mat = bfMatSub(firstProduct, secondProduct);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  bfMatDelete(&firstProduct);
  bfMatDelete(&secondProduct);

  return mat;
}

BfMat *bfMatDiffToType(BfMatDiff const *matDiff, BfType type) {
  BF_ERROR_BEGIN();

  BfMat *firstConverted = bfMatToType(matDiff->first, type);
  HANDLE_ERROR();

  BfMat *secondConverted = bfMatToType(matDiff->second, type);
  HANDLE_ERROR();

  BfMat *result = bfMatSub(firstConverted, secondConverted);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  bfMatDelete(&firstConverted);
  bfMatDelete(&secondConverted);

  return result;
}

/** Upcasting: MatDiff -> Mat */

BfMat *bfMatDiffToMat(BfMatDiff *matDiff) {
  return &matDiff->super;
}

BfMat const *bfMatDiffConstToMatConst(BfMatDiff const  *matDiff) {
  return &matDiff->super;
}

/** Downcasting: Mat -> MatDiff */

BfMatDiff const *bfMatConstToMatDiffConst(BfMat const *mat) {
  BF_ERROR_BEGIN();

  BfMatDiff const *matDiff = NULL;

  if (bfMatGetType(mat) != BF_TYPE_MAT_DIFF)
    RAISE_ERROR(BF_ERROR_TYPE_ERROR);

  matDiff = (BfMatDiff const *)mat;

  BF_ERROR_END() {
    BF_DIE();
  }

  return matDiff;
}

/** Implementation: MatDiff */

BfMatDiff *bfMatDiffAlloc(void) {
  BF_ERROR_BEGIN();

  BfMatDiff *matDiff = bfMemAlloc(1, sizeof(BfMatDiff));
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  return matDiff;
}

BfMatDiff *bfMatDiffNew(BfMat *first, BfMat *second, BfPolicy policy) {
  BF_ERROR_BEGIN();

  BfMatDiff *matDiff = bfMatDiffAlloc();
  HANDLE_ERROR();

  bfMatDiffInit(matDiff, first, second, policy);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  return matDiff;
}

void bfMatDiffInit(BfMatDiff *matDiff, BfMat *first, BfMat *second, BfPolicy policy) {
  BF_ERROR_BEGIN();

  BfSize numRows = bfMatGetNumRows(first);
  if (numRows != bfMatGetNumRows(second))
    RAISE_ERROR(BF_ERROR_INCOMPATIBLE_SHAPES);

  BfSize numCols = bfMatGetNumCols(first);
  if (numCols != bfMatGetNumCols(second))
    RAISE_ERROR(BF_ERROR_INCOMPATIBLE_SHAPES);

  bfMatInit(&matDiff->super, &MAT_VTABLE, numRows, numCols);
  HANDLE_ERROR();

  matDiff->first = bfMatGet(first, policy);
  HANDLE_ERROR();

  matDiff->second = bfMatGet(second, policy);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }
}

void bfMatDiffDeinit(BfMatDiff *matDiff) {
  bfMatDelete(&matDiff->first);
  bfMatDelete(&matDiff->second);
  bfMatDeinit(&matDiff->super);
}

void bfMatDiffDealloc(BfMatDiff **matDiff) {
  bfMemFree(*matDiff);
  *matDiff = NULL;
}

void bfMatDiffDeinitAndDealloc(BfMatDiff **matDiff) {
  bfMatDiffDeinit(*matDiff);
  bfMatDiffDealloc(matDiff);
}

BfMat *bfMatDiffGetFirst(BfMatDiff *matDiff) {
  return matDiff->first;
}

BfMat *bfMatDiffGetSecond(BfMatDiff *matDiff) {
  return matDiff->second;
}
