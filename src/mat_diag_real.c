#include <bf/mat_diag_real.h>

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <bf/error.h>
#include <bf/error_macros.h>

#define INTERFACE BF_INTERFACE_Mat
BF_DEFINE_VTABLE(Mat, MatDiagReal);
#undef INTERFACE

BfMat *bfMatDiagRealToMat(BfMatDiagReal *mat) {
  return &mat->super;
}

BfMat const *bfMatDiagRealConstToMatConst(BfMatDiagReal const *mat) {
  return &mat->super;
}

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

  bfMatDiagRealToMat(view)->props |= BF_MAT_PROPS_VIEW;

  END_ERROR_HANDLING() {
    free(view);
    view = NULL;
  }

  return view;
}

void bfMatDiagRealInit(BfMatDiagReal *mat, BfSize numRows, BfSize numCols) {
  BEGIN_ERROR_HANDLING();

  bfMatInit(&mat->super, &MatVtbl, numRows, numCols);
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

void bfMatDiagRealDealloc(BfMatDiagReal **mat) {
  free(*mat);
  *mat = NULL;
}

void bfMatDiagRealDeinitAndDealloc(BfMatDiagReal **mat) {
  bfMatDiagRealDeinit(*mat);
  bfMatDiagRealDealloc(mat);
}

BfMatDiagReal *
bfMatDiagRealGetDiagBlock(BfMatDiagReal *mat, BfSize i0, BfSize i1) {
  assert(i0 < i1);
  assert(i1 <= mat->numElts);

  BfMat *super = bfMatDiagRealToMat(mat);
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
    bfMatDiagRealDeinitAndDealloc(&submat);

  return submat;
}

BfMatDenseComplex *
bfMatDiagRealDenseComplexSolve(BfMatDiagReal const *lhs,
                               BfMatDenseComplex const *rhs)
{
  BEGIN_ERROR_HANDLING();

  BfMat const *lhsSuper = bfMatDiagRealConstToMatConst(lhs);
  BfMat const *rhsSuper = bfMatDenseComplexConstToMatConst(rhs);

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
    bfMatDenseComplexDeinitAndDealloc(&result);

  return result;
}

/* Interface: Mat */

void bfMatDiagRealDelete(BfMat **mat) {
  bfMatDiagRealDeinitAndDealloc((BfMatDiagReal **)mat);
}

BfMat *bfMatDiagRealEmptyLike(BfMat const *, BfSize, BfSize) {
  assert(false);
}

BfMat *bfMatDiagRealZerosLike(BfMat const *, BfSize, BfSize) {
  assert(false);
}

BfMatType bfMatDiagRealGetType(BfMat const *mat) {
  (void)mat;
  return BF_MAT_TYPE_DIAG_REAL;
}

bool bfMatDiagRealInstanceOf(BfMat const *mat, BfMatType matType) {
  return bfMatTypeDerivedFrom(bfMatGetType(mat), matType);
}

BfSize bfMatDiagRealNumBytes(BfMat const *mat) {
  (void)mat;
  assert(false);
  return BF_SIZE_BAD_VALUE;
}

void bfMatDiagRealSave(BfMat const *mat, char const *path) {
  (void)mat;
  (void)path;
  assert(false);
}

BfSize bfMatDiagRealGetNumRows(BfMat const *mat) { assert(false); }

BfSize bfMatDiagRealGetNumCols(BfMat const *mat) { assert(false); }

BfMat *bfMatDiagRealGetRowRange(BfMat *mat, BfSize i0, BfSize i1) {
  assert(false);
}

BfMat *bfMatDiagRealGetColRange(BfMat *mat, BfSize j0, BfSize j1) {
  assert(false);
}

void bfMatDiagRealSetRowRange(BfMat *mat, BfSize i0, BfSize i1,
                              BfMat const *otherMat) {
  assert(false);
}

void bfMatDiagRealAddInplace(BfMat *mat, BfMat const *otherMat) {
  assert(false);
}

BfMat *bfMatDiagRealMul(BfMat const *mat, BfMat const *otherMat) {
  assert(false);
}

BfMat *bfMatDiagRealLstSq(BfMat const *mat, BfMat const *otherMat) {
  assert(false);
}
