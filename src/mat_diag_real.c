#include <bf/mat_diag_real.h>

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mat_coo_real.h>

#define INTERFACE BF_INTERFACE_Mat
BF_DEFINE_VTABLE(Mat, MatDiagReal)
#undef INTERFACE

/* Interface: Mat */

BF_STUB(BfMat *, MatDiagRealCopy, BfMat const *)

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

BF_STUB(BfVec *, MatDiagRealGetRowCopy, BfMat const *, BfSize)
BF_STUB(BfVec *, MatDiagRealGetRowView, BfMat *, BfSize)
BF_STUB(BfVec *, MatDiagRealGetColView, BfMat *, BfSize)
BF_STUB(BfVec *, MatDiagRealGetColRangeView, BfMat *, BfSize, BfSize, BfSize)

void bfMatDiagRealDelete(BfMat **mat) {
  bfMatDiagRealDeinitAndDealloc((BfMatDiagReal **)mat);
}

BF_STUB(BfMat *, MatDiagRealEmptyLike, BfMat const *, BfSize, BfSize)
BF_STUB(BfMat *, MatDiagRealZerosLike, BfMat const *, BfSize, BfSize)

BfType bfMatDiagRealGetType(BfMat const *mat) {
  (void)mat;
  return BF_TYPE_MAT_DIAG_REAL;
}

bool bfMatDiagRealInstanceOf(BfMat const *mat, BfType type) {
  return bfTypeDerivedFrom(bfMatGetType(mat), type);
}

BF_STUB(BfSize, MatDiagRealNumBytes, BfMat const *)
BF_STUB(void, MatDiagRealSave, BfMat const *, char const *)
BF_STUB(void, MatDiagRealPrint, BfMat const *, FILE *)

BfSize bfMatDiagRealGetNumRows(BfMat const *mat) {
  if (!bfMatInstanceOf(mat, BF_TYPE_MAT_DIAG_REAL)) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return BF_SIZE_BAD_VALUE;
  } else {
    return mat->numRows;
  }
}

BfSize bfMatDiagRealGetNumCols(BfMat const *mat) {
  if (!bfMatInstanceOf(mat, BF_TYPE_MAT_DIAG_REAL)) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return BF_SIZE_BAD_VALUE;
  } else {
    return mat->numCols;
  }
}

BF_STUB(void, MatDiagRealSetRow, BfMat *, BfSize, BfVec const *)
BF_STUB(void, MatDiagRealSetCol, BfMat *, BfSize, BfVec const *)
BF_STUB(void, MatDiagRealSetColRange, BfMat *, BfSize, BfSize, BfSize, BfVec const *)
BF_STUB(BfMat *, MatDiagRealGetRowRange, BfMat *, BfSize, BfSize)
BF_STUB(BfMat *, MatDiagRealGetColRange, BfMat *, BfSize, BfSize)

BfMat *bfMatDiagRealGetRowRangeCopy(BfMat const *mat, BfSize i0, BfSize i1) {
  BEGIN_ERROR_HANDLING();

  BfMatDiagReal const *matDiagReal = NULL;
  BfMatCooReal *rowRange = NULL;

  if (i0 > i1)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  matDiagReal = bfMatConstToMatDiagRealConst(mat);
  HANDLE_ERROR();

  if (i1 > mat->numRows)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  rowRange = bfMatCooRealNew();
  HANDLE_ERROR();

  /* Need to make sure we get the correct number of elements in the
   * COO matrix when we slice */
  BfSize numElts = matDiagReal->numElts;
  if (i1 - i0 < numElts)
    numElts = i1 - i0;

  bfMatCooRealInitEmpty(rowRange, i1 - i0, mat->numCols, numElts);
  HANDLE_ERROR();

  /* Fill the COO matrix */
  BfReal *ptr = matDiagReal->data;
  for (BfSize k = 0; k < numElts; ++k) {
    rowRange->rowInd[k] = k;
    rowRange->colInd[k] = i0 + k;
    rowRange->value[k] = *ptr++;
  }

  END_ERROR_HANDLING()
    bfMatCooRealDeinitAndDealloc(&rowRange);

  return bfMatCooRealToMat(rowRange);
}

BF_STUB(BfMat *, MatDiagRealGetColRangeCopy, BfMat const *, BfSize, BfSize)
BF_STUB(void, MatDiagRealSetRowRange, BfMat *, BfSize, BfSize, BfMat const *)
BF_STUB(void, MatDiagRealPermuteRows, BfMat *, BfPerm const *)
BF_STUB(void, MatDiagRealPermuteCols, BfMat *, BfPerm const *)
BF_STUB(BfVec *, MatDiagRealRowDists, BfMat const *, BfMat const *)
BF_STUB(BfVec *, MatDiagRealColDists, BfMat const *, BfMat const *)
BF_STUB(BfVec *, MatDiagRealColDots, BfMat const *, BfMat const *)
BF_STUB(BfVec *, MatDiagRealColNorms, BfMat const *)
BF_STUB(void, MatDiagRealScaleCols, BfMat *, BfVec const *)
BF_STUB(BfVec *, MatDiagRealSumCols, BfMat const *)
BF_STUB(void, MatDiagRealAddInplace, BfMat *, BfMat const *)
BF_STUB(void, MatDiagRealAddDiag, BfMat *, BfMat const *)
BF_STUB(BfMat *, MatDiagRealSub, BfMat const *, BfMat const *)
BF_STUB(void, MatDiagRealSubInplace, BfMat *, BfMat const *)
BF_STUB(BfMat *, MatDiagRealMul, BfMat const *, BfMat const *)
BF_STUB(BfVec *, MatDiagRealMulVec, BfMat const *, BfVec const *)
BF_STUB(void, MatDiagRealMulInplace, BfMat *, BfMat const *)
BF_STUB(BfMat *, MatDiagRealSolveLU, BfMat const *, BfMat const *)
BF_STUB(BfMat *, MatDiagRealLstSq, BfMat const *, BfMat const *)
BF_STUB(BfMat *, MatDiagRealGetGivensRotation, BfVec const *, BfSize, BfSize)
BF_STUB(bool, MatDiagRealIsUpperTri, BfMat const *)
BF_STUB(BfVec *, MatDiagRealBackwardSolveVec, BfMat const *, BfVec const *)
BF_STUB(bool, MatDiagRealIsZero, BfMat const *)
BF_STUB(void, MatDiagRealNegate, BfMat *)
BF_STUB(BfMat *, MatDiagRealToType, BfMat const *, BfType)

/* Upcasting: */

BfMat *bfMatDiagRealToMat(BfMatDiagReal *mat) {
  return &mat->super;
}

BfMat const *bfMatDiagRealConstToMatConst(BfMatDiagReal const *mat) {
  return &mat->super;
}

/* Downcasting: */

BfMatDiagReal *bfMatToMatDiagReal(BfMat *mat) {
  if (!bfMatInstanceOf(mat, BF_TYPE_MAT_DIAG_REAL)) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return NULL;
  } else {
    return (BfMatDiagReal *)mat;
  }
}

BfMatDiagReal const *bfMatConstToMatDiagRealConst(BfMat const *mat) {
  if (!bfMatInstanceOf(mat, BF_TYPE_MAT_DIAG_REAL)) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return NULL;
  } else {
    return (BfMatDiagReal const *)mat;
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

void bfMatDiagRealInitView(BfMatDiagReal *mat, BfSize numRows, BfSize numCols,
                           BfReal *data) {
  BEGIN_ERROR_HANDLING();

  bfMatInit(&mat->super, &MatVtbl, numRows, numCols);
  HANDLE_ERROR();

  mat->super.props |= BF_MAT_PROPS_VIEW;

  mat->numElts = numRows < numCols ? numRows : numCols;
  mat->data = data;

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

  BfMatDiagReal *matView = bfMatDiagRealNew();
  HANDLE_ERROR();

  BfSize numElts = i1 - i0;
  BfReal *data = mat->data + i0;

  bfMatDiagRealInitView(matView, numElts, numElts, data);

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
