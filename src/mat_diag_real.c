#include <bf/mat_diag_real.h>

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <bf/error.h>
#include <bf/error_macros.h>

#define INTERFACE BF_INTERFACE_Mat
BF_DEFINE_VTABLE(Mat, MatDiagReal)
#undef INTERFACE

/* Interface: Mat */

BfMat *bfMatDiagRealGetView(BfMat *mat) {
  BEGIN_ERROR_HANDLING();

  BfMatDiagReal *view = malloc(sizeof(BfMatDiagReal));
  if (view == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BfMatDiagReal *matDiagReal = bfMatToMatDiagReal(mat);

  *view = *matDiagReal;

  bfMatDiagRealToMat(view)->props |= BF_MAT_PROPS_VIEW;

  END_ERROR_HANDLING() {
    free(view);
    view = NULL;
  }

  return bfMatDiagRealToMat(view);
}

void bfMatDiagRealDelete(BfMat **mat) {
  bfMatDiagRealDeinitAndDealloc((BfMatDiagReal **)mat);
}

BF_STUB(BfMat *, MatDiagRealEmptyLike, BfMat const *, BfSize, BfSize)
BF_STUB(BfMat *, MatDiagRealZerosLike, BfMat const *, BfSize, BfSize)

BfMatType bfMatDiagRealGetType(BfMat const *mat) {
  (void)mat;
  return BF_MAT_TYPE_DIAG_REAL;
}

bool bfMatDiagRealInstanceOf(BfMat const *mat, BfMatType matType) {
  return bfMatTypeDerivedFrom(bfMatGetType(mat), matType);
}

BF_STUB(BfSize, MatDiagRealNumBytes, BfMat const *)
BF_STUB(void, MatDiagRealSave, BfMat const *, char const *)
BF_STUB(void, MatDiagRealPrint, FILE *, BfMat const *)
BF_STUB(BfSize, MatDiagRealGetNumRows, BfMat const *)
BF_STUB(BfSize, MatDiagRealGetNumCols, BfMat const *)
BF_STUB(BfMat *, MatDiagRealGetRowRange, BfMat *, BfSize, BfSize)
BF_STUB(BfMat *, MatDiagRealGetColRange, BfMat *, BfSize, BfSize)
BF_STUB(void, MatDiagRealSetRowRange, BfMat *, BfSize, BfSize, BfMat const *)
BF_STUB(BfMat *, MatDiagRealRowDists, BfMat const *, BfMat const *)
BF_STUB(void, MatDiagRealScaleCols, BfMat *, BfMat const *)
BF_STUB(BfMat *, MatDiagRealSumCols, BfMat const *)
BF_STUB(void, MatDiagRealAddInplace, BfMat *, BfMat const *)
BF_STUB(void, MatDiagRealAddDiag, BfMat *, BfMat const *)
BF_STUB(BfMat *, MatDiagRealMul, BfMat const *, BfMat const *)
BF_STUB(void, MatDiagRealMulInplace, BfMat *, BfMat const *)
BF_STUB(BfMat *, MatDiagRealSolve, BfMat const *, BfMat const *)
BF_STUB(BfMat *, MatDiagRealLstSq, BfMat const *, BfMat const *)

/* Upcasting: */

BfMat *bfMatDiagRealToMat(BfMatDiagReal *mat) {
  return &mat->super;
}

BfMat const *bfMatDiagRealConstToMatConst(BfMatDiagReal const *mat) {
  return &mat->super;
}

/* Downcasting: */

BfMatDiagReal *bfMatToMatDiagReal(BfMat *mat) {
  if (!bfMatInstanceOf(mat, BF_MAT_TYPE_DIAG_REAL)) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return NULL;
  } else {
    return (BfMatDiagReal *)mat;
  }
}

/* Implementation: MatDiagReal */

BfMatDiagReal *bfMatDiagRealNew() {
  BEGIN_ERROR_HANDLING();

  BfMatDiagReal *mat = malloc(sizeof(BfMatDiagReal));
  if (mat == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING() {}

  return mat;
}

BfMatDiagReal *bfMatDiagRealEye(BfSize numRows, BfSize numCols) {
  BEGIN_ERROR_HANDLING();

  BfMatDiagReal *mat = bfMatDiagRealNew();
  HANDLE_ERROR();

  bfMatDiagRealInit(mat, numRows, numCols);
  HANDLE_ERROR();

  bfMatDiagRealSetConstant(mat, 1);

  END_ERROR_HANDLING()
    mat = NULL;

  return mat;
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

void bfMatDiagRealSetConstant(BfMatDiagReal *mat, BfReal value) {
  for (BfSize i = 0; i < mat->numElts; ++i)
    mat->data[i] = value;
}

BfMatDiagReal *
bfMatDiagRealGetDiagBlock(BfMatDiagReal *mat, BfSize i0, BfSize i1) {
  assert(i0 < i1);
  assert(i1 <= mat->numElts);

  BfMat *super = bfMatDiagRealToMat(mat);
  assert(!bfMatIsTransposed(super)); // TODO: implement

  BEGIN_ERROR_HANDLING();

  BfMatDiagReal *matView = (BfMatDiagReal *)mat;
  HANDLE_ERROR();

  matView->numElts = i1 - i0;

  if (matView->numElts == mat->numElts)
    return matView;

  matView->super.numRows = matView->numElts;
  matView->super.numCols = matView->numElts;

  matView->data += i0;

  END_ERROR_HANDLING()
    bfMatDiagRealDeinitAndDealloc(&matView);

  return matView;
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
