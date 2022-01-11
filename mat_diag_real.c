#include "mat_diag_real.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "error.h"
#include "error_macros.h"

static BfMatVtable matDiagRealVtbl = {
  .deinit = (__typeof__(&bfMatDeinit))bfMatDiagRealDeinit,
  .delete = (__typeof__(&bfMatDelete))bfMatDiagRealDelete,
  .deinitAndDelete = (__typeof__(&bfMatDeinitAndDelete))bfMatDiagRealDeinitAndDelete,
  .getType = (__typeof__(&bfMatGetType))bfMatDiagRealGetType,
  .numBytes = (__typeof__(&bfMatNumBytes))bfMatDiagRealNumBytes,
  .save = (__typeof__(&bfMatSave))bfMatDiagRealSave,
  .mul = NULL, // (__typeof__(&bfMatMul))bfMatDiagRealMul,
  .lstSq = NULL, // (__typeof__(&bfMatLstSq))bfMatDiagRealLstSq
};

BfMatDiagReal *bfMatDiagRealNew() {
  BEGIN_ERROR_HANDLING();

  BfMatDiagReal *mat = malloc(sizeof(BfMatDiagReal));
  if (mat == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING() {}

  return mat;
}

BfMatDiagReal *bfMatDiagRealNewView(BfMatDiagReal *mat) {
  BEGIN_ERROR_HANDLING();

  BfMatDiagReal *view = malloc(sizeof(BfMatDiagReal));
  if (view == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  *view = *mat;

  bfMatDiagRealGetMatPtr(view)->props |= BF_MAT_PROPS_VIEW;

  END_ERROR_HANDLING() {
    free(view);
    view = NULL;
  }

  return view;
}

void bfMatDiagRealInit(BfMatDiagReal *mat, BfSize numRows, BfSize numCols) {
  BEGIN_ERROR_HANDLING();

  bfMatInit(&mat->super, &matDiagRealVtbl, numRows, numCols);
  HANDLE_ERROR();

  mat->numElts = numRows < numCols ? numRows : numCols;

  mat->data = malloc(mat->numElts*sizeof(BfReal));
  if (mat->data == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING()
    bfMatDiagRealDeinit(mat);
}

void bfMatDiagRealDeinit(BfMatDiagReal *mat) {
  mat->numElts = BF_SIZE_BAD_VALUE;

  if (!(mat->super.props & BF_MAT_PROPS_VIEW))
    free(mat->data);

  mat->data = NULL;
}

void bfMatDiagRealDelete(BfMatDiagReal **mat) {
  free(*mat);
  *mat = NULL;
}

void bfMatDiagRealDeinitAndDelete(BfMatDiagReal **mat) {
  bfMatDiagRealDeinit(*mat);
  bfMatDiagRealDelete(mat);
}

BfMat *bfMatDiagRealGetMatPtr(BfMatDiagReal *mat) {
  return &mat->super;
}

BfMat const *bfMatDiagRealGetMatConstPtr(BfMatDiagReal const *mat) {
  return &mat->super;
}

BfMatType bfMatDiagRealGetType(BfMatDiagReal const *mat) {
  (void)mat;
  return BF_MAT_TYPE_DIAG_REAL;
}

BfSize bfMatDiagRealNumBytes(BfMatDiagReal const *mat) {
  (void)mat;
  assert(false);
  return BF_SIZE_BAD_VALUE;
}

void bfMatDiagRealSave(BfMatDiagReal const *mat, char const *path) {
  (void)mat;
  (void)path;
  assert(false);
}

BfMatDiagReal *
bfMatDiagRealGetDiagBlock(BfMatDiagReal *mat, BfSize i0, BfSize i1) {
  assert(i0 < i1);
  assert(i1 <= mat->numElts);

  BfMat *super = bfMatDiagRealGetMatPtr(mat);
  assert(!bfMatIsTransposed(super)); // TODO: implement

  BEGIN_ERROR_HANDLING();

  BfMatDiagReal *submat = bfMatDiagRealNewView(mat);
  HANDLE_ERROR();

  submat->numElts = i1 - i0;

  if (submat->numElts == mat->numElts)
    return submat;

  submat->super.numRows = submat->numElts;
  submat->super.numCols = submat->numElts;

  submat->data += i0;

  END_ERROR_HANDLING()
    bfMatDiagRealDeinitAndDelete(&submat);

  return submat;
}

BfMatDenseComplex *
bfMatDiagRealDenseComplexSolve(BfMatDiagReal const *lhs,
                               BfMatDenseComplex const *rhs)
{
  BEGIN_ERROR_HANDLING();

  BfMat const *lhsSuper = bfMatDiagRealGetMatConstPtr(lhs);
  BfMat const *rhsSuper = bfMatDenseComplexGetMatConstPtr(rhs);

  BfSize k = rhsSuper->numRows;
  if (k != lhsSuper->numRows)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize m = lhsSuper->numCols;
  BfSize n = rhsSuper->numCols;

  BfMatDenseComplex *result = bfMatDenseComplexNew();
  if (result == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  bfMatDenseComplexInit(result, m, n);
  HANDLE_ERROR();

  BfReal const *lhsPtr = lhs->data;
  BfComplex const *rhsPtr = rhs->data;
  BfComplex *resultPtr = result->data;

  for (BfSize i = 0; i < k; ++i) {
    for (BfSize j = 0; j < n; ++j)
      *resultPtr++ = *rhsPtr++ / *lhsPtr;
    lhsPtr++;
  }

  memset(resultPtr, 0x0, (m - k)*n*sizeof(BfComplex));

  END_ERROR_HANDLING()
    bfMatDenseComplexDeinitAndDelete(&result);

  return result;
}
