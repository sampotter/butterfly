#include <bf/mat_diag_real.h>

#include <bf/assert.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mem.h>
#include <bf/mat_coo_real.h>
#include <bf/vec_real.h>

/** Interface: Mat */

static BfMatVtable MAT_VTABLE = {
  .GetView = (__typeof__(&bfMatDiagRealGetView))bfMatDiagRealGetView,
  .GetRowCopy = (__typeof__(&bfMatDiagRealGetRowCopy))bfMatDiagRealGetRowCopy,
  .Delete = (__typeof__(&bfMatDiagRealDelete))bfMatDiagRealDelete,
  .GetType = (__typeof__(&bfMatDiagRealGetType))bfMatDiagRealGetType,
  .GetNumRows = (__typeof__(&bfMatGetNumRows))bfMatDiagRealGetNumRows,
  .GetNumCols = (__typeof__(&bfMatGetNumCols))bfMatDiagRealGetNumCols,
  .GetRowRangeCopy = (__typeof__(&bfMatDiagRealGetRowRangeCopy))bfMatDiagRealGetRowRangeCopy,
  .MulVec = (__typeof__(&bfMatMulVec))bfMatDiagRealMulVec,
};

BfMat *bfMatDiagRealGetView(BfMat *mat) {
  BF_ERROR_BEGIN();

  BfMatDiagReal *view = bfMemAlloc(1, sizeof(BfMatDiagReal));
  if (view == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BfMatDiagReal *matDiagReal = bfMatToMatDiagReal(mat);

  *view = *matDiagReal;

  bfMatDiagRealToMat(view)->props |= BF_MAT_PROPS_VIEW;

  BF_ERROR_END() {
    bfMemFree(view);
    view = NULL;
  }

  return bfMatDiagRealToMat(view);
}

BfVec *bfMatDiagRealGetRowCopy(BfMat const *mat, BfSize i) {
  BF_ERROR_BEGIN();

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

  BF_ERROR_END()
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

BfSize bfMatDiagRealGetNumRows(BfMatDiagReal const *matDiagReal) {
  return matDiagReal->super.numRows;
}

BfSize bfMatDiagRealGetNumCols(BfMatDiagReal const *matDiagReal) {
  return matDiagReal->super.numCols;
}

BfMat *bfMatDiagRealGetRowRangeCopy(BfMat const *mat, BfSize i0, BfSize i1) {
  BF_ERROR_BEGIN();

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

  BF_ERROR_END()
    bfMatCooRealDeinitAndDealloc(&rowRange);

  return bfMatCooRealToMat(rowRange);
}

BfVec *mulVec_vecReal(BfMatDiagReal const *matDiagReal, BfVecReal const *vecReal) {
  BF_ERROR_BEGIN();

  BfVecReal *result = NULL;

  BfSize m = bfMatDiagRealGetNumRows(matDiagReal);
  BfSize n = bfMatDiagRealGetNumCols(matDiagReal);

  if (vecReal->super.size != n)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  result = bfVecRealNewWithValue(m, 0);

  BfReal *writePtr = result->data;
  BfReal const *scalePtr = matDiagReal->data;
  BfReal const *readPtr = vecReal->data;

  for (BfSize j = 0; j < n; ++j) {
    *writePtr = *readPtr;
    *writePtr *= *scalePtr;
    writePtr += result->stride;
    ++scalePtr;
    readPtr += vecReal->stride;
  }

  BF_ERROR_END() {
    BF_DIE();
  }

  return bfVecRealToVec(result);
}

BfVec *bfMatDiagRealMulVec(BfMatDiagReal const *matDiagReal, BfVec const *vec) {
  switch (bfVecGetType(vec)) {
  case BF_TYPE_VEC_REAL:
    return mulVec_vecReal(matDiagReal, bfVecConstToVecRealConst(vec));
  default:
    bfSetError(BF_ERROR_NOT_IMPLEMENTED);
    return NULL;
  }
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
  BF_ERROR_BEGIN();

  BfMatDiagReal *mat = bfMemAlloc(1, sizeof(BfMatDiagReal));
  if (mat == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BF_ERROR_END() {}

  return mat;
}

BfMatDiagReal *bfMatDiagRealEye(BfSize numRows, BfSize numCols) {
  return bfMatDiagRealNewConstant(numRows, numCols, 1.0);
}

BfMatDiagReal *bfMatDiagRealNewConstant(BfSize numRows, BfSize numCols, BfReal diagValue) {
  BF_ERROR_BEGIN();

  BfMatDiagReal *mat = bfMatDiagRealNew();
  HANDLE_ERROR();

  bfMatDiagRealInit(mat, numRows, numCols);
  HANDLE_ERROR();

  bfMatDiagRealSetConstant(mat, diagValue);

  BF_ERROR_END()
    mat = NULL;

  return mat;
}

BfMatDiagReal *bfMatDiagRealNewFromData(BfSize numRows, BfSize numCols, BfReal const *data) {
  BF_ERROR_BEGIN();

  BfMatDiagReal *matDiagReal = bfMatDiagRealNew();
  HANDLE_ERROR();

  bfMatDiagRealInit(matDiagReal, numRows, numCols);
  HANDLE_ERROR();

  bfMemCopy(data, matDiagReal->numElts, sizeof(BfReal), matDiagReal->data);

  BF_ERROR_END() {
    BF_DIE();
  }

  return matDiagReal;
}

void bfMatDiagRealInit(BfMatDiagReal *mat, BfSize numRows, BfSize numCols) {
  BF_ERROR_BEGIN();

  bfMatInit(&mat->super, &MAT_VTABLE, numRows, numCols);
  HANDLE_ERROR();

  mat->numElts = numRows < numCols ? numRows : numCols;

  mat->data = bfMemAlloc(mat->numElts, sizeof(BfReal));
  if (mat->data == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BF_ERROR_END()
    bfMatDiagRealDeinit(mat);
}

void bfMatDiagRealInitView(BfMatDiagReal *mat, BfSize numRows, BfSize numCols,
                           BfReal *data) {
  BF_ERROR_BEGIN();

  bfMatInit(&mat->super, &MAT_VTABLE, numRows, numCols);
  HANDLE_ERROR();

  mat->super.props |= BF_MAT_PROPS_VIEW;

  mat->numElts = numRows < numCols ? numRows : numCols;
  mat->data = data;

  BF_ERROR_END()
    bfMatDiagRealDeinit(mat);
}

void bfMatDiagRealDeinit(BfMatDiagReal *mat) {
  mat->numElts = BF_SIZE_BAD_VALUE;

  if (!(mat->super.props & BF_MAT_PROPS_VIEW))
    bfMemFree(mat->data);

  mat->data = NULL;
}

void bfMatDiagRealDealloc(BfMatDiagReal **mat) {
  bfMemFree(*mat);
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
  BF_ASSERT(i0 < i1);
  BF_ASSERT(i1 <= mat->numElts);

  BfMat *super = bfMatDiagRealToMat(mat);
  BF_ASSERT(!bfMatIsTransposed(super)); // TODO: implement

  BF_ERROR_BEGIN();

  BfMatDiagReal *matView = bfMatDiagRealNew();
  HANDLE_ERROR();

  BfSize numElts = i1 - i0;
  BfReal *data = mat->data + i0;

  bfMatDiagRealInitView(matView, numElts, numElts, data);

  BF_ERROR_END()
    bfMatDiagRealDeinitAndDealloc(&matView);

  return matView;
}

BfMatDenseComplex *
bfMatDiagRealDenseComplexSolve(BfMatDiagReal const *lhs,
                               BfMatDenseComplex const *rhs)
{
  BF_ERROR_BEGIN();

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

  bfMemZero(resultPtr, (m - k)*n, sizeof(BfComplex));

  BF_ERROR_END()
    bfMatDenseComplexDeinitAndDealloc(&result);

  return result;
}

BfVec *bfMatDiagRealGetVecView(BfMatDiagReal *mat) {
  BF_ERROR_BEGIN();

  BfVecReal *vecReal = bfVecRealNew();
  HANDLE_ERROR();

  bfVecRealInitView(vecReal, mat->numElts, 1, mat->data);

  BF_ERROR_END()
    bfVecRealDeinitAndDealloc(&vecReal);

  return bfVecRealToVec(vecReal);
}
