#include <bf/mat_diag_real.h>

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mat_coo_real.h>
#include <bf/vec_real.h>

/** Interface: Mat */

static BfMatVtable MAT_VTABLE = {
  .GetView = (__typeof__(&bfMatDiagRealGetView))bfMatDiagRealGetView,
  .GetRowCopy = (__typeof__(&bfMatDiagRealGetRowCopy))bfMatDiagRealGetRowCopy,
  .Delete = (__typeof__(&bfMatDiagRealDelete))bfMatDiagRealDelete,
  .GetType = (__typeof__(&bfMatDiagRealGetType))bfMatDiagRealGetType,
  .GetNumRows = (__typeof__(&bfMatDiagRealGetNumRows))bfMatDiagRealGetNumRows,
  .GetNumCols = (__typeof__(&bfMatDiagRealGetNumCols))bfMatDiagRealGetNumCols,
  .GetRowRangeCopy = (__typeof__(&bfMatDiagRealGetRowRangeCopy))bfMatDiagRealGetRowRangeCopy,
};

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

BfVec *bfMatDiagRealGetRowCopy(BfMat const *mat, BfSize i) {
  BEGIN_ERROR_HANDLING();

  BfSize numRows = bfMatGetNumRows(mat);
  if (i >= numRows)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize numCols = bfMatGetNumCols(mat);

  BfMatDiagReal const *matDiagReal = bfMatConstToMatDiagRealConst(mat);
  HANDLE_ERROR();

  BfVecReal *rowCopy = bfVecRealNew();
  HANDLE_ERROR();

  bfVecRealInit(rowCopy, numCols);
  HANDLE_ERROR();

  /* Fill row with zeros */
  BfReal *ptr = rowCopy->data;
  for (BfSize j = 0; j < numCols; ++j) {
    *ptr = 0;
    ptr += rowCopy->stride;
  }

  /* Set ith element to corresponding diagonal entry of mat */
  *(rowCopy->data + i*rowCopy->stride) = *(matDiagReal->data + i);

  END_ERROR_HANDLING()
    bfVecRealDeinitAndDealloc(&rowCopy);

  return bfVecRealToVec(rowCopy);
}

void bfMatDiagRealDelete(BfMat **mat) {
  bfMatDiagRealDeinitAndDealloc((BfMatDiagReal **)mat);
}

BfType bfMatDiagRealGetType(BfMat const *mat) {
  (void)mat;
  return BF_TYPE_MAT_DIAG_REAL;
}

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
  return bfMatDiagRealNewConstant(numRows, numCols, 1.0);
}

BfMatDiagReal *bfMatDiagRealNewConstant(BfSize numRows, BfSize numCols, BfReal diagValue) {
  BEGIN_ERROR_HANDLING();

  BfMatDiagReal *mat = bfMatDiagRealNew();
  HANDLE_ERROR();

  bfMatDiagRealInit(mat, numRows, numCols);
  HANDLE_ERROR();

  bfMatDiagRealSetConstant(mat, diagValue);

  END_ERROR_HANDLING()
    mat = NULL;

  return mat;
}

void bfMatDiagRealInit(BfMatDiagReal *mat, BfSize numRows, BfSize numCols) {
  BEGIN_ERROR_HANDLING();

  bfMatInit(&mat->super, &MAT_VTABLE, numRows, numCols);
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

  bfMatInit(&mat->super, &MAT_VTABLE, numRows, numCols);
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

BfVec *bfMatDiagRealGetVecView(BfMatDiagReal *mat) {
  BEGIN_ERROR_HANDLING();

  BfVecReal *vecReal = bfVecRealNew();
  HANDLE_ERROR();

  bfVecRealInitView(vecReal, mat->numElts, 1, mat->data);

  END_ERROR_HANDLING()
    bfVecRealDeinitAndDealloc(&vecReal);

  return bfVecRealToVec(vecReal);
}
