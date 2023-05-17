#include <bf/mat_givens.h>

#include <math.h>
#include <stdlib.h>

#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mem.h>

/** Interface: Mat */

static BfMatVtable MAT_VTABLE = {
  .GetType = (__typeof__(&bfMatGivensComplexGetType))bfMatGivensComplexGetType,
  .GetNumRows = (__typeof__(&bfMatGivensComplexGetNumRows))bfMatGivensComplexGetNumRows,
  .GetNumCols = (__typeof__(&bfMatGivensComplexGetNumCols))bfMatGivensComplexGetNumCols,
};

BfType bfMatGivensComplexGetType(BfMat const *mat) {
  (void)mat;
  return BF_TYPE_MAT_GIVENS_COMPLEX;
}

BfSize bfMatGivensComplexGetNumRows(BfMat const *mat) {
  if (!bfMatInstanceOf(mat, BF_TYPE_MAT_GIVENS_COMPLEX)) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return BF_SIZE_BAD_VALUE;
  } else {
    return mat->numRows;
  }
}

BfSize bfMatGivensComplexGetNumCols(BfMat const *mat) {
  if (!bfMatInstanceOf(mat, BF_TYPE_MAT_GIVENS_COMPLEX)) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return BF_SIZE_BAD_VALUE;
  } else {
    return mat->numCols;
  }
}

/** Upcasting: */

BfMat *bfMatGivensComplexToMat(BfMatGivensComplex *mat) {
  return &mat->super;
}

/** Downcasting: */

BfMatGivensComplex const *bfMatConstToMatGivensComplexConst(BfMat const *mat) {
  if (!bfMatInstanceOf(mat, BF_TYPE_MAT_GIVENS_COMPLEX)) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return NULL;
  } else {
    return (BfMatGivensComplex const *)mat;
  }
}

/** Implementation: MatGivensComplex */

BfMatGivensComplex *bfMatGivensComplexNew() {
  BF_ERROR_BEGIN();

  BfMatGivensComplex *mat = bfMemAlloc(1, sizeof(BfMatGivensComplex));
  if (mat == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING() {}

  return mat;
}

void bfMatGivensComplexInit(BfMatGivensComplex *mat, BfSize n,
                            BfSize srcInd, BfSize elimInd,
                            BfComplex c, BfComplex s) {
  BF_ERROR_BEGIN();

  bfMatInit(&mat->super, &MAT_VTABLE, n, n);
  HANDLE_ERROR();

  mat->srcInd = srcInd;
  mat->elimInd = elimInd;
  mat->c = c;
  mat->s = s;

  END_ERROR_HANDLING()
    bfMatGivensComplexDeinit(mat);
}

void bfMatGivensComplexDeinit(BfMatGivensComplex *mat) {
  mat->srcInd = BF_SIZE_BAD_VALUE;
  mat->elimInd = BF_SIZE_BAD_VALUE;
  mat->c = NAN + I*NAN;
  mat->s = NAN + I*NAN;
}

void bfMatGivensComplexDealloc(BfMatGivensComplex **mat) {
  free(*mat);
  *mat = NULL;
}

void bfMatGivensComplexDeinitAndDealloc(BfMatGivensComplex **mat) {
  bfMatGivensComplexDeinit(*mat);
  bfMatGivensComplexDealloc(mat);
}
