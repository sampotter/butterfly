#include <bf/mat_dense_complex.h>

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <openblas/lapacke.h>

#include <bf/blas.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/lu_dense_complex.h>
#include <bf/mat_coo_complex.h>
#include <bf/mat_coo_real.h>
#include <bf/mat_dense_real.h>
#include <bf/mat_diag_real.h>
#include <bf/vec_complex.h>
#include <bf/vec_real.h>

/** Static functions: */

static enum CBLAS_TRANSPOSE getCblasTranspose(BfMatDenseComplex const *mat) {
  BfMat const *super = bfMatDenseComplexConstToMatConst(mat);
  if (super->props & (BF_MAT_PROPS_TRANS | BF_MAT_PROPS_CONJ))
    return CblasConjTrans;
  else if (super->props & BF_MAT_PROPS_TRANS)
    return CblasTrans;
  else
    return CblasNoTrans;
}

static BfSize getLeadingDimension(BfMatDenseComplex const *mat) {
  return mat->rowStride;
}

static BfVec *bfMatDenseComplexDenseComplexColDists(
  BfMatDenseComplex const *matDenseComplex,
  BfMatDenseComplex const *otherMatDenseComplex)
{
  BEGIN_ERROR_HANDLING();

  BfVecReal *result = NULL;

  BfMat const *mat = bfMatDenseComplexConstToMatConst(matDenseComplex);
  BfMat const *otherMat = bfMatDenseComplexConstToMatConst(otherMatDenseComplex);

  BfSize numRows = bfMatDenseComplexGetNumRows(mat);
  if (numRows != bfMatDenseComplexGetNumRows(otherMat))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize numCols = bfMatDenseComplexGetNumCols(mat);
  if (numCols != bfMatDenseComplexGetNumCols(otherMat))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  result = bfVecRealNew();
  HANDLE_ERROR();

  bfVecRealInit(result, numCols);
  HANDLE_ERROR();

  BfComplex const *colPtr, *otherColPtr;
  BfReal *resultPtr = result->data;

  for (BfSize j = 0; j < numCols; ++j) {
    colPtr = matDenseComplex->data + j*matDenseComplex->colStride;
    otherColPtr = otherMatDenseComplex->data + j*otherMatDenseComplex->colStride;

    BfReal sum = 0;
    for (BfSize i = 0; i < numRows; ++i) {
      BfComplex diff = *otherColPtr - *colPtr;
      sum += diff*conj(diff);
      colPtr += matDenseComplex->rowStride;
      otherColPtr += otherMatDenseComplex->rowStride;
    }
    sum = sqrt(sum);

    *resultPtr++ = sum;
  }

  END_ERROR_HANDLING()
    result = NULL;

  return bfVecRealToVec(result);
}

static BfVec *bfMatDenseComplexDenseComplexColDots(
  BfMatDenseComplex const *matDenseComplex,
  BfMatDenseComplex const *otherMatDenseComplex)
{
  BEGIN_ERROR_HANDLING();

  BfVecComplex *result = NULL;

  BfMat const *mat = bfMatDenseComplexConstToMatConst(matDenseComplex);
  BfMat const *otherMat = bfMatDenseComplexConstToMatConst(otherMatDenseComplex);

  BfSize numRows = bfMatDenseComplexGetNumRows(mat);
  if (numRows != bfMatDenseComplexGetNumRows(otherMat))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize numCols = bfMatDenseComplexGetNumCols(mat);
  if (numCols != bfMatDenseComplexGetNumCols(otherMat))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  result = bfVecComplexNew();
  HANDLE_ERROR();

  bfVecComplexInit(result, numCols);
  HANDLE_ERROR();

  BfComplex const *colPtr, *otherColPtr;
  BfComplex *resultPtr = result->data;

  for (BfSize j = 0; j < numCols; ++j) {
    colPtr = matDenseComplex->data + j*matDenseComplex->colStride;
    otherColPtr = otherMatDenseComplex->data + j*otherMatDenseComplex->colStride;

    BfComplex dot = cblas_zdotc(
      /* n: */ numRows,
      /* x: */ colPtr,
      /* incx: */ matDenseComplex->rowStride,
      /* y: */ otherColPtr,
      /* incy: */ otherMatDenseComplex->rowStride);

    *resultPtr++ = dot;
  }

  END_ERROR_HANDLING()
    result = NULL;

  return bfVecComplexToVec(result);
}

void bfMatDenseComplexScaleCols_real(BfMatDenseComplex *mat, BfVec const *vec) {
  BEGIN_ERROR_HANDLING();

  BfVecReal const *vecReal = bfVecConstToVecRealConst(vec);
  HANDLE_ERROR();

  for (BfSize i = 0; i < mat->super.numRows; ++i) {
    BfReal *vecRealData = vecReal->data;
    BfComplex *rowData = mat->data + i*mat->rowStride;
    for (BfSize j = 0; j < mat->super.numCols; ++j) {
      *(rowData + j*mat->colStride) *= *vecRealData;
      vecRealData += vecReal->stride;
    }
  }

  END_ERROR_HANDLING() {}
}

void bfMatDenseComplexScaleCols_complex(BfMatDenseComplex *mat, BfVec const *vec) {
  BEGIN_ERROR_HANDLING();

  BfVecComplex const *vecComplex = bfVecConstToVecComplexConst(vec);
  HANDLE_ERROR();

  BfSize m = mat->super.numRows;
  BfSize n = mat->super.numCols;

  for (BfSize j = 0; j < n; ++j) {
    BfComplex const *readPtr = vecComplex->data + j*vecComplex->stride;
    BfComplex *writePtr = mat->data + j*mat->colStride;
    cblas_zscal(m, readPtr, writePtr, mat->rowStride);
  }

  END_ERROR_HANDLING() {}
}

/** Interface: Mat */

static BfMatVtable MAT_VTABLE = {
  .Copy = (__typeof__(&bfMatDenseComplexCopy))bfMatDenseComplexCopy,
  .GetView = (__typeof__(&bfMatDenseComplexGetView))bfMatDenseComplexGetView,
  .GetRowCopy = (__typeof__(&bfMatDenseComplexGetRowCopy))bfMatDenseComplexGetRowCopy,
  .GetRowView = (__typeof__(&bfMatDenseComplexGetRowView))bfMatDenseComplexGetRowView,
  .GetColView = (__typeof__(&bfMatDenseComplexGetColView))bfMatDenseComplexGetColView,
  .GetColRangeView = (__typeof__(&bfMatDenseComplexGetColRangeView))bfMatDenseComplexGetColRangeView,
  .Delete = (__typeof__(&bfMatDenseComplexDelete))bfMatDenseComplexDelete,
  .EmptyLike = (__typeof__(&bfMatDenseComplexEmptyLike))bfMatDenseComplexEmptyLike,
  .ZerosLike = (__typeof__(&bfMatDenseComplexZerosLike))bfMatDenseComplexZerosLike,
  .GetType = (__typeof__(&bfMatDenseComplexGetType))bfMatDenseComplexGetType,
  .NumBytes = (__typeof__(&bfMatDenseComplexNumBytes))bfMatDenseComplexNumBytes,
  .Save = (__typeof__(&bfMatDenseComplexSave))bfMatDenseComplexSave,
  .Print = (__typeof__(&bfMatDenseComplexPrint))bfMatDenseComplexPrint,
  .GetNumRows = (__typeof__(&bfMatDenseComplexGetNumRows))bfMatDenseComplexGetNumRows,
  .GetNumCols = (__typeof__(&bfMatDenseComplexGetNumCols))bfMatDenseComplexGetNumCols,
  .SetRow = (__typeof__(&bfMatDenseComplexSetRow))bfMatDenseComplexSetRow,
  .SetCol = (__typeof__(&bfMatDenseComplexSetCol))bfMatDenseComplexSetCol,
  .SetColRange = (__typeof__(&bfMatDenseComplexSetColRange))bfMatDenseComplexSetColRange,
  .GetRowRange = (__typeof__(&bfMatDenseComplexGetRowRange))bfMatDenseComplexGetRowRange,
  .GetColRange = (__typeof__(&bfMatDenseComplexGetColRange))bfMatDenseComplexGetColRange,
  .SetRowRange = (__typeof__(&bfMatDenseComplexSetRowRange))bfMatDenseComplexSetRowRange,
  .PermuteRows = (__typeof__(&bfMatDenseComplexPermuteRows))bfMatDenseComplexPermuteRows,
  .ColDists = (__typeof__(&bfMatDenseComplexColDists))bfMatDenseComplexColDists,
  .ColDots = (__typeof__(&bfMatDenseComplexColDots))bfMatDenseComplexColDots,
  .ColNorms = (__typeof__(&bfMatDenseComplexColNorms))bfMatDenseComplexColNorms,
  .Scale = (__typeof__(&bfMatScale))bfMatDenseComplexScale,
  .ScaleCols = (__typeof__(&bfMatDenseComplexScaleCols))bfMatDenseComplexScaleCols,
  .AddInplace = (__typeof__(&bfMatDenseComplexAddInplace))bfMatDenseComplexAddInplace,
  .Sub = (__typeof__(&bfMatDenseComplexSub))bfMatDenseComplexSub,
  .SubInplace = (__typeof__(&bfMatDenseComplexSubInplace))bfMatDenseComplexSubInplace,
  .Mul = (__typeof__(&bfMatDenseComplexMul))bfMatDenseComplexMul,
  .MulVec = (__typeof__(&bfMatDenseComplexMulVec))bfMatDenseComplexMulVec,
  .Solve = (__typeof__(&bfMatSolve))bfMatDenseComplexSolve,
  .SolveLU = (__typeof__(&bfMatDenseComplexSolveLU))bfMatDenseComplexSolveLU,
  .LstSq = (__typeof__(&bfMatDenseComplexLstSq))bfMatDenseComplexLstSq,
  .IsUpperTri = (__typeof__(&bfMatDenseComplexIsUpperTri))bfMatDenseComplexIsUpperTri,
  .BackwardSolveVec = (__typeof__(&bfMatDenseComplexBackwardSolveVec))bfMatDenseComplexBackwardSolveVec,
  .Negate = (__typeof__(&bfMatDenseComplexNegate))bfMatDenseComplexNegate,
  .GetBlockView = (__typeof__(&bfMatGetBlockView))bfMatDenseComplexGetBlockView,
  .GetLu = (__typeof__(&bfMatGetLu))bfMatDenseComplexGetLu,
  .DivideCols = (__typeof__(&bfMatDivideCols))bfMatDenseComplexDivideCols,
};

BfMat *bfMatDenseComplexCopy(BfMat const *mat) {
  BEGIN_ERROR_HANDLING();

  BfMatDenseComplex const *matDenseComplex = NULL;
  BfMatDenseComplex *copy = NULL;

  matDenseComplex = bfMatConstToMatDenseComplexConst(mat);
  HANDLE_ERROR();

  BfSize m = bfMatGetNumRows(mat);
  BfSize n = bfMatGetNumCols(mat);

  copy = bfMatDenseComplexNew();
  HANDLE_ERROR();

  bfMatDenseComplexInit(copy, m, n);
  HANDLE_ERROR();

  for (BfSize i = 0; i < m; ++i) {
    BfComplex *outPtr = copy->data + i*copy->rowStride;
    BfComplex const *inPtr = matDenseComplex->data + i*matDenseComplex->rowStride;
    for (BfSize j = 0; j < n; ++j) {
      *outPtr = *inPtr;
      outPtr += copy->colStride;
      inPtr += matDenseComplex->colStride;
    }
  }

  END_ERROR_HANDLING()
    bfMatDenseComplexDeinitAndDealloc(&copy);

  return bfMatDenseComplexToMat(copy);
}

BfMat *bfMatDenseComplexGetView(BfMat *mat) {
  BEGIN_ERROR_HANDLING();

  BfMatDenseComplex *matDenseComplex = bfMatToMatDenseComplex(mat);
  HANDLE_ERROR();

  BfMatDenseComplex *matDenseComplexView = bfMatDenseComplexNew();

  *matDenseComplexView = *matDenseComplex;

  BfMat *matView = bfMatDenseComplexToMat(matDenseComplexView);

  matView->props |= BF_MAT_PROPS_VIEW;

  END_ERROR_HANDLING()
    matView = NULL;

  return matView;
}

BfVec *bfMatDenseComplexGetRowCopy(BfMat const *mat, BfSize i) {
  BEGIN_ERROR_HANDLING();

  BfMatDenseComplex const *matDenseComplex = NULL;
  BfVecComplex *rowCopy = NULL;

  matDenseComplex = bfMatConstToMatDenseComplexConst(mat);
  HANDLE_ERROR();

  BfSize numRows = bfMatGetNumRows(mat);
  if (i >= numRows)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize numCols = bfMatGetNumCols(mat);

  rowCopy = bfVecComplexNew();
  HANDLE_ERROR();

  bfVecComplexInit(rowCopy, numCols);
  HANDLE_ERROR();

  BfComplex *writePtr = rowCopy->data;
  BfComplex const *readPtr = matDenseComplex->data + i*matDenseComplex->rowStride;
  for (BfSize j = 0; j < numCols; ++j) {
    *writePtr = *readPtr;
    writePtr += rowCopy->stride;
    readPtr += matDenseComplex->colStride;
  }

  END_ERROR_HANDLING() {}

  return bfVecComplexToVec(rowCopy);
}

BfVec *bfMatDenseComplexGetRowView(BfMat *mat, BfSize i) {
  BEGIN_ERROR_HANDLING();

  BfMatDenseComplex *matDenseComplex = bfMatToMatDenseComplex(mat);
  HANDLE_ERROR();

  BfVecComplex *rowView = bfVecComplexNew();
  HANDLE_ERROR();

  BfSize n = bfMatGetNumCols(mat);
  BfSize stride = matDenseComplex->rowStride;
  BfComplex *data = matDenseComplex->data + i*matDenseComplex->rowStride;

  bfVecComplexInitView(rowView, n, stride, data);

  END_ERROR_HANDLING()
    bfVecComplexDeinitAndDealloc(&rowView);

  return bfVecComplexToVec(rowView);
}

BfVec *bfMatDenseComplexGetColView(BfMat *mat, BfSize j) {
  BEGIN_ERROR_HANDLING();

  BfMatDenseComplex *matDenseComplex = bfMatToMatDenseComplex(mat);
  HANDLE_ERROR();

  BfVecComplex *colView = bfVecComplexNew();
  HANDLE_ERROR();

  BfSize m = bfMatGetNumRows(mat);
  BfSize stride = matDenseComplex->rowStride;
  BfComplex *data = matDenseComplex->data + j*matDenseComplex->colStride;

  bfVecComplexInitView(colView, m, stride, data);

  END_ERROR_HANDLING()
    bfVecComplexDeinitAndDealloc(&colView);

  return bfVecComplexToVec(colView);
}

BfVec *bfMatDenseComplexGetColRangeView(BfMat *mat, BfSize i0, BfSize i1, BfSize j) {
  BEGIN_ERROR_HANDLING();

  BfMatDenseComplex *matDenseComplex = bfMatToMatDenseComplex(mat);
  HANDLE_ERROR();

  BfVecComplex *colView = bfVecComplexNew();
  HANDLE_ERROR();

  if (i0 > i1)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize m = bfMatGetNumRows(mat);
  if (i1 > m)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize stride = matDenseComplex->rowStride;
  BfComplex *data = matDenseComplex->data + j*matDenseComplex->colStride
    + i0*matDenseComplex->rowStride;

  bfVecComplexInitView(colView, i1 - i0, stride, data);

  END_ERROR_HANDLING()
    bfVecComplexDeinitAndDealloc(&colView);

  return bfVecComplexToVec(colView);
}

void bfMatDenseComplexDelete(BfMat **mat) {
  bfMatDenseComplexDeinitAndDealloc((BfMatDenseComplex **)mat);
}

BfMat *bfMatDenseComplexEmptyLike(BfMat const *mat, BfSize numRows, BfSize numCols) {
  BEGIN_ERROR_HANDLING();

  BfMatDenseComplex *result = NULL;

  if (numRows == BF_SIZE_BAD_VALUE)
    numRows = bfMatDenseComplexGetNumRows(mat);

  if (numCols == BF_SIZE_BAD_VALUE)
    numCols = bfMatDenseComplexGetNumCols(mat);

  result = bfMatDenseComplexNew();
  HANDLE_ERROR();

  bfMatDenseComplexInit(result, numRows, numCols);
  HANDLE_ERROR();

  END_ERROR_HANDLING()
    bfMatDenseComplexDeinitAndDealloc(&result);

  return bfMatDenseComplexToMat(result);
}

BfMat *bfMatDenseComplexZerosLike(BfMat const *mat, BfSize numRows, BfSize numCols) {
  BEGIN_ERROR_HANDLING();

  BfMatDenseComplex *zeros = NULL;

  if (numRows == BF_SIZE_BAD_VALUE)
    numRows = bfMatDenseComplexGetNumRows(mat);

  if (numCols == BF_SIZE_BAD_VALUE)
    numCols = bfMatDenseComplexGetNumCols(mat);

  zeros = bfMatDenseComplexZeros(numRows, numCols);
  HANDLE_ERROR();

  END_ERROR_HANDLING()
    bfMatDenseComplexDeinitAndDealloc(&zeros);

  return bfMatDenseComplexToMat(zeros);
}

BfType bfMatDenseComplexGetType(BfMat const *mat) {
  (void)mat;
  return BF_TYPE_MAT_DENSE_COMPLEX;
}

BfSize bfMatDenseComplexNumBytes(BfMat const *mat) {
  (void)mat;
  assert(false);
  return BF_SIZE_BAD_VALUE;
}

void bfMatDenseComplexSave(BfMat const *mat, char const *path) {
  BEGIN_ERROR_HANDLING();

  BfMatDenseComplex const *matDenseComplex = bfMatConstToMatDenseComplexConst(mat);

  FILE *fp = fopen(path, "w");
  if (fp == NULL)
    RAISE_ERROR(BF_ERROR_FILE_ERROR);

  for (BfSize i = 0; i < mat->numRows; ++i) {
    BfComplex const *ptr = matDenseComplex->data + i*matDenseComplex->rowStride;
    for (BfSize j = 0; j < mat->numCols; ++j) {
      fwrite(ptr, sizeof(BfComplex), 1, fp);
      ptr += matDenseComplex->colStride;
    }
  }

  END_ERROR_HANDLING() {
    assert(false);
  }

  fclose(fp);
}

void bfMatDenseComplexPrint(BfMat const *mat, FILE *fp) {
  BEGIN_ERROR_HANDLING();

  BfMatDenseComplex const *matDenseComplex = bfMatConstToMatDenseComplexConst(mat);
  HANDLE_ERROR();

  BfSize m = bfMatGetNumRows(mat);
  BfSize n = bfMatGetNumCols(mat);

  for (BfSize i = 0; i < m; ++i) {
    BfComplex const *rowPtr = matDenseComplex->data + i*matDenseComplex->rowStride;
    for (BfSize j = 0; j < n; ++j) {
      BfComplex z = *rowPtr;
      fprintf(fp, " %g + i*%g", creal(z), cimag(z));
      rowPtr += matDenseComplex->colStride;
    }
    fprintf(fp, "\n");
  }

  END_ERROR_HANDLING() {}
}

BfSize bfMatDenseComplexGetNumRows(BfMat const *mat) {
  return bfMatIsTransposed(mat) ? mat->numCols : mat->numRows;
}

BfSize bfMatDenseComplexGetNumCols(BfMat const *mat) {
  return bfMatIsTransposed(mat) ? mat->numRows : mat->numCols;
}

static void
setRow_complex(BfMatDenseComplex *matDenseComplex, BfSize i,
               BfVecComplex const *rowVec, BfSize n) {
  BfComplex *outPtr = matDenseComplex->data + i*matDenseComplex->rowStride;
  BfComplex const *inPtr = rowVec->data;
  for (BfSize j = 0; j < n; ++j) {
    *outPtr = *inPtr;
    outPtr += matDenseComplex->colStride;
    inPtr += rowVec->stride;
  }
}

static void
setRow_real(BfMatDenseComplex *matDenseComplex, BfSize i,
            BfVecReal const *rowVec, BfSize n) {
  BfComplex *outPtr = matDenseComplex->data + i*matDenseComplex->rowStride;
  BfReal const *inPtr = rowVec->data;
  for (BfSize j = 0; j < n; ++j) {
    *outPtr = *inPtr;
    outPtr += matDenseComplex->colStride;
    inPtr += rowVec->stride;
  }
}

void bfMatDenseComplexSetRow(BfMat *mat, BfSize i, BfVec const *rowVec) {
  BEGIN_ERROR_HANDLING();

  BfMatDenseComplex *matDenseComplex = bfMatToMatDenseComplex(mat);
  HANDLE_ERROR();

  if (i >= bfMatGetNumRows(mat))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize n = bfMatGetNumCols(mat);
  if (rowVec->size != n)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  switch (bfVecGetType(rowVec)) {
  case BF_TYPE_VEC_COMPLEX:
    setRow_complex(matDenseComplex, i, bfVecConstToVecComplexConst(rowVec), n);
    break;
  case BF_TYPE_VEC_REAL:
    setRow_real(matDenseComplex, i, bfVecConstToVecRealConst(rowVec), n);
    break;
  default:
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);
  }

  END_ERROR_HANDLING() {}
}

void bfMatDenseComplexSetCol(BfMat *mat, BfSize j, BfVec const *vec) {
  BEGIN_ERROR_HANDLING();

  BfMatDenseComplex *matDenseComplex = NULL;
  BfVecComplex const *vecComplex = NULL;

  BfSize m = bfMatDenseComplexGetNumRows(mat);
  if (m != vec->size)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  matDenseComplex = bfMatToMatDenseComplex(mat);
  HANDLE_ERROR();

  vecComplex = bfVecConstToVecComplexConst(vec);
  HANDLE_ERROR();

  BfComplex *outPtr = matDenseComplex->data + j*matDenseComplex->colStride;
  BfComplex const *inPtr = vecComplex->data;
  for (BfSize i = 0; i < m; ++i) {
    *outPtr = *inPtr;
    outPtr += matDenseComplex->rowStride;
    inPtr += vecComplex->stride;
  }

  END_ERROR_HANDLING() {}
}

void bfMatDenseComplexSetColRange(BfMat *mat, BfSize j, BfSize i0, BfSize i1,
                                  BfVec const *vec) {
  BEGIN_ERROR_HANDLING();

  BfMatDenseComplex *matDenseComplex = NULL;
  BfVecComplex const *vecComplex = NULL;

  if (i0 > i1)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize m = bfMatDenseComplexGetNumRows(mat);
  if (m < i1 - i0)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  matDenseComplex = bfMatToMatDenseComplex(mat);
  HANDLE_ERROR();

  vecComplex = bfVecConstToVecComplexConst(vec);
  HANDLE_ERROR();

  BfComplex *outPtr = matDenseComplex->data
    + i0*matDenseComplex->rowStride + j*matDenseComplex->colStride;
  BfComplex const *inPtr = vecComplex->data + i0*vecComplex->stride;
  for (BfSize i = i0; i < i1; ++i) {
    *outPtr = *inPtr;
    outPtr += matDenseComplex->rowStride;
    inPtr += vecComplex->stride;
  }

  END_ERROR_HANDLING() {}
}

BfMat *bfMatDenseComplexGetRowRange(BfMat *mat, BfSize i0, BfSize i1) {
  BfSize numRows = mat->numRows;

  assert(i0 < i1);
  assert(i1 <= numRows);
  assert(!bfMatIsTransposed(mat)); // TODO: implement

  BEGIN_ERROR_HANDLING();

  BfMat *matView = bfMatGetView(mat);

  BfMatDenseComplex *submat = bfMatToMatDenseComplex(matView);
  HANDLE_ERROR();

  if (i1 - i0 != numRows) {
    bfMatDenseComplexToMat(submat)->numRows = i1 - i0;
    submat->data += submat->rowStride*i0;
  }

  END_ERROR_HANDLING()
    bfMatDelete(&matView);

  return matView;
}

BfMat *bfMatDenseComplexGetColRange(BfMat *mat, BfSize j0, BfSize j1) {
  BEGIN_ERROR_HANDLING();

  if (bfMatIsTransposed(mat))
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);

  BfSize numCols = mat->numCols;
  if (!(j0 <= j1 && j1 <= numCols))
    RAISE_ERROR(BF_ERROR_OUT_OF_RANGE);

  BfMat *matView = bfMatGetView(mat);

  BfMatDenseComplex *matDenseComplexView = bfMatToMatDenseComplex(matView);
  HANDLE_ERROR();

  if (j1 - j0 < numCols) {
    matDenseComplexView->super.numCols = j1 - j0;
    matDenseComplexView->data += matDenseComplexView->colStride*j0;
  }

  END_ERROR_HANDLING()
    bfMatDenseComplexDeinitAndDealloc(&matDenseComplexView);

  return bfMatDenseComplexToMat(matDenseComplexView);
}

void bfMatDenseComplexSetRowRange(BfMat *mat, BfSize i0, BfSize i1, BfMat const *rows) {
  BEGIN_ERROR_HANDLING();

  BfMatDenseComplex *matRows =
    bfMatToMatDenseComplex(bfMatDenseComplexGetRowRange(mat, i0, i1));
  HANDLE_ERROR();

  switch (bfMatGetType(rows)) {
  case BF_TYPE_MAT_DENSE_COMPLEX:
    bfMatDenseComplexSet(matRows, bfMatConstToMatDenseComplexConst(rows));
    HANDLE_ERROR();
    break;
  default:
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);
  }

  END_ERROR_HANDLING() {}
}

void bfMatDenseComplexPermuteRows(BfMat *mat, BfPerm const *perm) {
  BEGIN_ERROR_HANDLING();

  BfSize m = bfMatGetNumRows(mat);
  BfSize n = bfMatGetNumCols(mat);

  BfMatDenseComplex *matDenseComplex = NULL;
  BfMatDenseComplex *matDenseComplexPerm = NULL;

  BfMat *matPerm = bfMatCopy(mat);
  HANDLE_ERROR();

  matDenseComplex = bfMatToMatDenseComplex(mat);
  HANDLE_ERROR();

  matDenseComplexPerm = bfMatToMatDenseComplex(matPerm);
  HANDLE_ERROR();

  for (BfSize i = 0; i < m; ++i) {
    BfComplex const *inRowPtr =
      matDenseComplexPerm->data + i*matDenseComplexPerm->rowStride;
    BfComplex *outRowPtr =
      matDenseComplex->data + perm->index[i]*matDenseComplex->rowStride;
    for (BfSize j = 0; j < n; ++j) {
      *outRowPtr = *inRowPtr;
      inRowPtr += matDenseComplexPerm->colStride;
      outRowPtr += matDenseComplex->colStride;
    }
  }

  END_ERROR_HANDLING() {}

  bfMatDelete(&matPerm);
}

BfVec *bfMatDenseComplexColDists(BfMat const *mat, BfMat const *otherMat) {
  BEGIN_ERROR_HANDLING();

  BfMatDenseComplex const *matDenseComplex = bfMatConstToMatDenseComplexConst(mat);
  HANDLE_ERROR();

  BfVec *result = NULL;

  switch (bfMatGetType(otherMat)) {
  case BF_TYPE_MAT_DENSE_COMPLEX:
    result = bfMatDenseComplexDenseComplexColDists(
        matDenseComplex, bfMatConstToMatDenseComplexConst(otherMat));
    break;
  default:
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);
  }

  END_ERROR_HANDLING()
    result = NULL;

  return result;
}

BfVec *bfMatDenseComplexColDots(BfMat const *mat, BfMat const *otherMat) {
  BEGIN_ERROR_HANDLING();

  BfMatDenseComplex const *matDenseComplex = bfMatConstToMatDenseComplexConst(mat);
  HANDLE_ERROR();

  BfVec *result = NULL;

  switch (bfMatGetType(otherMat)) {
  case BF_TYPE_MAT_DENSE_COMPLEX:
    result = bfMatDenseComplexDenseComplexColDots(
        matDenseComplex, bfMatConstToMatDenseComplexConst(otherMat));
    break;
  default:
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);
  }

  END_ERROR_HANDLING()
    result = NULL;

  return result;
}

BfVec *bfMatDenseComplexColNorms(BfMat const *mat) {
  BEGIN_ERROR_HANDLING();

  BfMatDenseComplex const *matDenseComplex = bfMatConstToMatDenseComplexConst(mat);
  HANDLE_ERROR();

  BfSize m = bfMatGetNumRows(mat);
  BfSize n = bfMatGetNumCols(mat);

  BfVecReal *colNorms = bfVecRealNew();
  HANDLE_ERROR();

  bfVecRealInit(colNorms, n);

  for (BfSize j = 0; j < n; ++j) {
    BfComplex const *readPtr = matDenseComplex->data + j*matDenseComplex->colStride;
    BfReal *writePtr = colNorms->data + j*colNorms->stride;
    *writePtr = cblas_dznrm2(m, readPtr, matDenseComplex->rowStride);
  }

  END_ERROR_HANDLING()
    bfVecRealDeinitAndDealloc(&colNorms);

  return bfVecRealToVec(colNorms);
}

void bfMatDenseComplexScale(BfMatDenseComplex *matDenseComplex, BfComplex scalar) {
  BfMat *mat = bfMatDenseComplexToMat(matDenseComplex);

  BfSize m = bfMatGetNumRows(mat);
  BfSize n = bfMatGetNumCols(mat);

  for (BfSize i = 0; i < m; ++i) {
    BfComplex *rowPtr = matDenseComplex->data + i*matDenseComplex->rowStride;
    for (BfSize j = 0; j < n; ++j) {
      *rowPtr *= scalar;
      rowPtr += matDenseComplex->colStride;
    }
  }
}

void bfMatDenseComplexScaleCols(BfMat *mat, BfVec const *vec) {
  BEGIN_ERROR_HANDLING();

  if (bfMatGetNumCols(mat) != vec->size)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfMatDenseComplex *matDenseComplex = bfMatToMatDenseComplex(mat);
  HANDLE_ERROR();

  BfType type = bfVecGetType(vec);

  switch (type) {
  case BF_TYPE_VEC_REAL:
    bfMatDenseComplexScaleCols_real(matDenseComplex, vec);
    break;
  case BF_TYPE_VEC_COMPLEX:
    bfMatDenseComplexScaleCols_complex(matDenseComplex, vec);
    break;
  default:
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);
  }

  END_ERROR_HANDLING() {}
}

void bfMatDenseComplexAddInplace(BfMat *mat, BfMat const *otherMat) {
  BEGIN_ERROR_HANDLING();

  BfMatDenseComplex *matDenseComplex = bfMatToMatDenseComplex(mat);

  switch (bfMatGetType(otherMat)) {
  case BF_TYPE_MAT_COO_COMPLEX:
    bfMatDenseComplexCooComplexAddInplace(
      matDenseComplex, bfMatConstToMatCooComplexConst(otherMat));
    break;
  case BF_TYPE_MAT_COO_REAL:
    bfMatDenseComplexCooRealAddInplace(
      matDenseComplex, bfMatConstToMatCooRealConst(otherMat));
    break;
  case BF_TYPE_MAT_DENSE_COMPLEX:
    bfMatDenseComplexDenseComplexAddInplace(
      matDenseComplex, bfMatConstToMatDenseComplexConst(otherMat));
    break;
  case BF_TYPE_MAT_DIAG_REAL:
    bfMatDenseComplexDiagRealAddInplace(
      matDenseComplex, bfMatConstToMatDiagRealConst(otherMat));
    break;
  default:
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);
  }

  END_ERROR_HANDLING() {}
}

BfMat *bfMatDenseComplexSub(BfMat const *mat, BfMat const *otherMat) {
  BEGIN_ERROR_HANDLING();

  BfMat *result = NULL;

  BfMatDenseComplex const *matDenseComplex = bfMatConstToMatDenseComplexConst(mat);
  (void)matDenseComplex;

  switch (bfMatGetType(otherMat)) {
  case BF_TYPE_MAT_DENSE_COMPLEX:
    result = bfMatDenseComplexToMat(
      bfMatDenseComplexDenseComplexSub(
        matDenseComplex, bfMatConstToMatDenseComplexConst(otherMat)));
    break;
  default:
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);
  }

  END_ERROR_HANDLING() {}

  return result;
}

static void
subInplace_denseComplex(BfMatDenseComplex *matDenseComplex,
                        BfMatDenseComplex const *otherMatDenseComplex,
                        BfSize m, BfSize n) {
  for (BfSize i = 0; i < m; ++i) {
    BfComplex *outPtr = matDenseComplex->data + i*matDenseComplex->rowStride;
    BfComplex const *inPtr = otherMatDenseComplex->data + i*otherMatDenseComplex->rowStride;
    for (BfSize j = 0; j < n; ++j) {
      *outPtr -= *inPtr;
      outPtr += matDenseComplex->colStride;
      inPtr += otherMatDenseComplex->colStride;
    }
  }
}

void bfMatDenseComplexSubInplace(BfMat *mat, BfMat const *otherMat) {
  BEGIN_ERROR_HANDLING();

  BfMatDenseComplex *matDenseComplex = bfMatToMatDenseComplex(mat);
  HANDLE_ERROR();

  BfSize m = bfMatGetNumRows(mat);
  if (m != bfMatGetNumRows(otherMat))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize n = bfMatGetNumCols(mat);
  if (n != bfMatGetNumCols(otherMat))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  switch (bfMatGetType(otherMat)) {
  case BF_TYPE_MAT_DENSE_COMPLEX:
    subInplace_denseComplex(
      matDenseComplex, bfMatConstToMatDenseComplexConst(otherMat), m, n);
    break;
  default:
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);
  }

  END_ERROR_HANDLING() {}
}

BfMat *bfMatDenseComplexMul(BfMat const *op1, BfMat const *op2) {
  BEGIN_ERROR_HANDLING();

  BfMatDenseComplex const *matDenseComplex = bfMatConstToMatDenseComplexConst(op1);

  BfMat *result = NULL;

  switch (bfMatGetType(op2)) {
  case BF_TYPE_MAT_DENSE_COMPLEX:
    result = bfMatDenseComplexToMat(
      bfMatDenseComplexDenseComplexMul(
        matDenseComplex, bfMatConstToMatDenseComplexConst(op2)));
    break;
  default:
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);
  }

  END_ERROR_HANDLING()
    bfMatDelete(&result);

  return result;
}

static BfVecComplex *
mulVec_complex(BfMatDenseComplex const *matDenseComplex,
               BfVecComplex const *vecComplex,
               BfSize m, BfSize n)
{
  BEGIN_ERROR_HANDLING();

  enum CBLAS_TRANSPOSE trans = getCblasTranspose(matDenseComplex);

  BfVecComplex *result = bfVecComplexNew();
  HANDLE_ERROR();

  bfVecComplexInit(result, m);
  HANDLE_ERROR();

  BfSize lda = getLeadingDimension(matDenseComplex);

  BfComplex alpha = 1, beta = 0;

  cblas_zgemv(CblasRowMajor, trans, m, n, &alpha, matDenseComplex->data,
              lda, vecComplex->data, vecComplex->stride, &beta, result->data,
              result->stride);

  END_ERROR_HANDLING()
    bfVecComplexDeinitAndDealloc(&result);

  return result;
}

BfVec *bfMatDenseComplexMulVec(BfMat const *mat, BfVec const *vec) {
  BEGIN_ERROR_HANDLING();

  BfMatDenseComplex const *matDenseComplex = NULL;
  BfVec *result = NULL;

  BfSize m = bfMatGetNumRows(mat);
  BfSize n = bfMatGetNumCols(mat);
  if (n != vec->size)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  matDenseComplex = bfMatConstToMatDenseComplexConst(mat);
  HANDLE_ERROR();

  switch (bfVecGetType(vec)) {
  case BF_TYPE_VEC_COMPLEX:
    result = bfVecComplexToVec(
      mulVec_complex(matDenseComplex, bfVecConstToVecComplexConst(vec), m, n));
    break;
  default:
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);
  }

  END_ERROR_HANDLING()
    bfVecDelete(&result);

  return result;
}

static BfMat *
solve_matDenseComplex_tri(BfMatDenseComplex const *matDenseComplex, BfMat const *otherMat) {
  BEGIN_ERROR_HANDLING();

  BfMat const *mat = bfMatDenseComplexConstToMatConst(matDenseComplex);

  BfSize m = bfMatGetNumRows(mat);
  BfSize n = bfMatGetNumCols(mat);
  BfSize p = bfMatGetNumCols(otherMat);

  if (m != n)
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);

  if (bfMatGetNumRows(otherMat) != m);

  assert(mat->props & BF_MAT_PROPS_TRI);

  if (mat->props & BF_MAT_PROPS_TRANS)
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);

  BfMatDenseComplex const *otherMatDenseComplex = bfMatConstToMatDenseComplexConst(otherMat);

  BfMatDenseComplex *resultMatDenseComplex = bfMatDenseComplexNew();
  HANDLE_ERROR();

  bfMatDenseComplexInit(resultMatDenseComplex, n, p);
  HANDLE_ERROR();

  memcpy(resultMatDenseComplex->data, otherMatDenseComplex->data, n*p*sizeof(BfComplex));

  BfComplex alpha = 1;

  cblas_ztrsm(
    CblasRowMajor,
    CblasLeft,
    mat->props & BF_MAT_PROPS_LOWER_TRI ? CblasLower : CblasUpper,
    CblasNoTrans,
    mat->props & BF_MAT_PROPS_UNIT ? CblasUnit : CblasNonUnit,
    bfMatGetNumRows(otherMat),
    bfMatGetNumCols(otherMat),
    &alpha,
    /* a: */ matDenseComplex->data,
    /* lda: */ n,
    /* b: */ resultMatDenseComplex->data,
    /* ldb: */ p);

  END_ERROR_HANDLING() {
    assert(false);
  }

  return bfMatDenseComplexToMat(resultMatDenseComplex);
}

static BfMat *solve_matDenseComplex(BfMatDenseComplex const *matDenseComplex, BfMat const *otherMat) {
  BfMat const *mat = bfMatDenseComplexConstToMatConst(matDenseComplex);

  if (mat->props & BF_MAT_PROPS_TRI)
    return solve_matDenseComplex_tri(matDenseComplex, otherMat);

  bfSetError(BF_ERROR_NOT_IMPLEMENTED);

  return NULL;
}

BfMat *bfMatDenseComplexSolve(BfMatDenseComplex const *matDenseComplex, BfMat const *otherMat) {
  switch (bfMatGetType(otherMat)) {
  case BF_TYPE_MAT_DENSE_COMPLEX:
    return solve_matDenseComplex(matDenseComplex, otherMat);
  default:
    bfSetError(BF_ERROR_NOT_IMPLEMENTED);
    return NULL;
  }
}

BfMat *bfMatDenseComplexSolveLU(BfMat const *A, BfMat const *B) {
  BEGIN_ERROR_HANDLING();

  BfMatDenseComplex const *matDenseComplexA = bfMatConstToMatDenseComplexConst(A);
  BfMat *result = NULL;

  switch (bfMatGetType(B)) {
  case BF_TYPE_MAT_DENSE_COMPLEX:
    result = bfMatDenseComplexToMat(
      bfMatDenseComplexDenseComplexSolve(
        matDenseComplexA, bfMatConstToMatDenseComplexConst(B)));
    break;
  default:
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);
  }

  END_ERROR_HANDLING()
    bfMatDelete(&result);

  return result;
}

BfMat *bfMatDenseComplexLstSq(BfMat const *lhs, BfMat const *rhs) {
  BEGIN_ERROR_HANDLING();

  BfMatDenseComplex const *matDenseComplex = bfMatConstToMatDenseComplexConst(lhs);

  BfMat *result;

  switch (bfMatGetType(rhs)) {
  case BF_TYPE_MAT_DENSE_COMPLEX:
    result = bfMatDenseComplexToMat(
      bfMatDenseComplexDenseComplexLstSq(
        matDenseComplex, bfMatConstToMatDenseComplexConst(rhs)));
    break;
  default:
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);
  }

  END_ERROR_HANDLING()
    bfMatDelete(&result);

  return result;
}

bool bfMatDenseComplexIsUpperTri(BfMat const *mat) {
  BEGIN_ERROR_HANDLING();

  BfMatDenseComplex const *matDenseComplex = NULL;

  matDenseComplex = bfMatConstToMatDenseComplexConst(mat);
  HANDLE_ERROR();

  BfSize m = bfMatGetNumRows(mat);

  bool upperTri = true;
  for (BfSize i = 0; i < m; ++i) {
    BfComplex const *rowPtr = matDenseComplex->data + i*matDenseComplex->rowStride;
    for (BfSize j = 0; j < i; ++j) {
      if (*rowPtr != 0) {
        upperTri = false;
        break;
      }
      rowPtr += matDenseComplex->colStride;
    }
  }

  END_ERROR_HANDLING() {}

  return upperTri;
}

static BfVec *
backwardSolveVec_complex(BfMatDenseComplex const *matDenseComplex,
                         BfVecComplex const *vecComplex, BfSize m)
{
  BEGIN_ERROR_HANDLING();

  BfVecComplex *result = bfVecToVecComplex(bfVecCopy(&vecComplex->super));
  HANDLE_ERROR();

  enum CBLAS_TRANSPOSE trans = getCblasTranspose(matDenseComplex);
  BfSize lda = getLeadingDimension(matDenseComplex);

  cblas_ztrsv(CblasRowMajor, CblasUpper, trans, CblasNonUnit, m,
              matDenseComplex->data, lda, result->data, result->stride);

  END_ERROR_HANDLING() {}

  return bfVecComplexToVec(result);
}

BfVec *bfMatDenseComplexBackwardSolveVec(BfMat const *mat, BfVec const *vec) {
  BEGIN_ERROR_HANDLING();

  BfMatDenseComplex const *matDenseComplex = NULL;
  BfVec *result = NULL;

  matDenseComplex = bfMatConstToMatDenseComplexConst(mat);
  HANDLE_ERROR();

  BfSize m = bfMatGetNumRows(mat);

  if (m != bfMatGetNumCols(mat))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (m != vec->size)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  switch (bfVecGetType(vec)) {
  case BF_TYPE_VEC_COMPLEX:
    result = backwardSolveVec_complex(
      matDenseComplex, bfVecConstToVecComplexConst(vec), m);
    HANDLE_ERROR();
    break;
  default:
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);
  }

  END_ERROR_HANDLING()
    bfVecDelete(&result);

  return result;
}

void bfMatDenseComplexNegate(BfMat *mat) {
  BEGIN_ERROR_HANDLING();

  BfMatDenseComplex *matDenseComplex = bfMatToMatDenseComplex(mat);
  HANDLE_ERROR();

  BfSize numRows = bfMatGetNumRows(mat);
  BfSize numCols = bfMatGetNumCols(mat);

  for (BfSize i = 0; i < numRows; ++i) {
    BfComplex *rowPtr = matDenseComplex->data + i*matDenseComplex->rowStride;
    for (BfSize j = 0; j < numCols; ++j) {
      *rowPtr *= -1;
      rowPtr += matDenseComplex->colStride;
    }
  }

  END_ERROR_HANDLING() {}
}

BfMat *bfMatDenseComplexGetBlockView(BfMatDenseComplex *mat, BfSize i0, BfSize i1, BfSize j0, BfSize j1) {
  BEGIN_ERROR_HANDLING();

  BfMatDenseComplex *blockView = bfMatDenseComplexNew();
  HANDLE_ERROR();

  bfMatInit(&blockView->super, &MAT_VTABLE, i1 - i0, j1 - j0);

  blockView->super.props |= BF_MAT_PROPS_VIEW;

  blockView->rowStride = mat->rowStride;
  blockView->colStride = mat->colStride;
  blockView->data = mat->data + i0*mat->rowStride + j0*mat->colStride;

  END_ERROR_HANDLING() {}

  return bfMatDenseComplexToMat(blockView);
}

BfLu *bfMatDenseComplexGetLu(BfMatDenseComplex const *matDenseComplex) {
  BEGIN_ERROR_HANDLING();

  BfMat const *mat = bfMatDenseComplexConstToMatConst(matDenseComplex);

  BfLuDenseComplex *luDenseComplex = bfLuDenseComplexNew();
  HANDLE_ERROR();

  bfLuDenseComplexInit(luDenseComplex, mat);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {
    assert(false);
  }

  return bfLuDenseComplexToLu(luDenseComplex);
}

static void divideCols_vecReal(BfMatDenseComplex *matDenseComplex, BfVec const *vec) {
  BEGIN_ERROR_HANDLING();

  BfMat *mat = bfMatDenseComplexToMat(matDenseComplex);

  BfSize numRows = bfMatGetNumRows(mat);
  BfSize numCols = bfMatGetNumCols(mat);

  if (numCols != vec->size)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfVecReal const *vecReal = bfVecConstToVecRealConst(vec);

  for (BfSize j = 0; j < numCols; ++j) {
    BfReal alpha = *(vecReal->data + j*vecReal->stride);
    BfComplex *writePtr = matDenseComplex->data + j*matDenseComplex->colStride;
    cblas_zdscal(numRows, 1/alpha, writePtr, matDenseComplex->rowStride);
  }

  END_ERROR_HANDLING() {
    assert(false);
  }
}

void bfMatDenseComplexDivideCols(BfMatDenseComplex *matDenseComplex, BfVec const *vec) {
  switch (bfVecGetType(vec)) {
  case BF_TYPE_VEC_REAL:
    divideCols_vecReal(matDenseComplex, vec);
    break;
  default:
    bfSetError(BF_ERROR_NOT_IMPLEMENTED);
  }
}

/* Implementation: MatDenseComplex */

BfSize bfMatDenseComplexGetRowStride(BfMatDenseComplex const *mat) {
  return bfMatIsTransposed(bfMatDenseComplexConstToMatConst(mat)) ?
    mat->colStride :
    mat->rowStride;
}

BfSize bfMatDenseComplexGetColStride(BfMatDenseComplex const *mat) {
  return bfMatIsTransposed(bfMatDenseComplexConstToMatConst(mat)) ?
    mat->rowStride :
    mat->colStride;
}

void bfMatDenseComplexSet(BfMatDenseComplex *dst, BfMatDenseComplex const *src) {
  BEGIN_ERROR_HANDLING();

  BfMat *mat = bfMatDenseComplexToMat(dst);
  BfMat const *matOther = bfMatDenseComplexConstToMatConst(src);

  BfSize dstRows = bfMatDenseComplexGetNumRows(mat);
  if (dstRows != bfMatDenseComplexGetNumRows(matOther))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize dstCols = bfMatDenseComplexGetNumCols(mat);
  if (dstCols != bfMatDenseComplexGetNumCols(matOther))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize dstRowStride = bfMatDenseComplexGetRowStride(dst);
  BfSize dstColStride = bfMatDenseComplexGetColStride(dst);

  BfSize srcRowStride = bfMatDenseComplexGetRowStride(src);
  BfSize srcColStride = bfMatDenseComplexGetColStride(src);

  for (BfSize i = 0; i < dstRows; ++i) {
    BfComplex *dstData = dst->data + i*dstRowStride;
    BfComplex const *srcData = src->data + i*srcRowStride;
    for (BfSize j = 0; j < dstCols; ++j)
      *(dstData + j*dstColStride) = *(srcData + j*srcColStride);
  }

  END_ERROR_HANDLING() {}
}

void bfMatDenseComplexSvd(BfMatDenseComplex const *mat, BfMatDenseComplex *U,
                          BfMatDiagReal *S, BfMatDenseComplex *VH) {
  BEGIN_ERROR_HANDLING();

  BfMat const *super = bfMatDenseComplexConstToMatConst(mat);

  BfSize m = super->numRows;
  BfSize n = super->numCols;

  BfReal *superb = NULL;

  /* zgesvd will overwrite A, so allocate space for a copy */
  BfComplex *dataCopy = malloc(m*n*sizeof(BfComplex));
  if (dataCopy == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  /* copy contents of A */
  memcpy(dataCopy, mat->data, m*n*sizeof(BfComplex));
  /* TODO: error-handling */

  /* output array which contains information about superdiagonal
   * elements which didn't converge
   *
   * more info here: tinyurl.com/2p8f5ev3 */
  superb = malloc((((m < n) ? m : n) - 1)*sizeof(BfReal));
  if (superb == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  // TODO: see call to dgesvd in mat_dense_real.c to figure out proper
  // way to set leading dimensions in call to zgesvd below...
  assert(m == n);

  /* compute the SVD */
  lapack_int info = LAPACKE_zgesvd(
    LAPACK_ROW_MAJOR, 'S', 'S', m, n, dataCopy, m, S->data, U->data, m,
    VH->data, n, superb);

  /* check for invalid arguments */
  if (info < 0)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  /* check for errors */
  if (info > 0)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  END_ERROR_HANDLING() {}

  bfMatDenseComplexToMat(U)->props |= BF_MAT_PROPS_ORTHO;
  bfMatDenseComplexToMat(VH)->props |= BF_MAT_PROPS_ORTHO;

  free(dataCopy);
  free(superb);
}

void bfMatDenseComplexDenseComplexAddInplace(BfMatDenseComplex *op1,
                                             BfMatDenseComplex const *op2)
{
  BEGIN_ERROR_HANDLING();

  BfMat *mat = bfMatDenseComplexToMat(op1);
  BfMat const *otherMat = bfMatDenseComplexConstToMatConst(op2);

  BfSize numRows = bfMatDenseComplexGetNumRows(mat);
  if (numRows != bfMatDenseComplexGetNumRows(otherMat))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize numCols = bfMatDenseComplexGetNumCols(mat);
  if (numCols != bfMatDenseComplexGetNumCols(otherMat))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  for (BfSize i = 0; i < numRows*numCols; ++i)
    op1->data[i] += op2->data[i];

  END_ERROR_HANDLING() {}
}

void bfMatDenseComplexDiagRealAddInplace(BfMatDenseComplex *op1,
                                         BfMatDiagReal const *op2) {
  BEGIN_ERROR_HANDLING();

  BfMat *mat = bfMatDenseComplexToMat(op1);
  BfMat const *otherMat = bfMatDiagRealConstToMatConst(op2);

  BfSize numRows = bfMatDenseComplexGetNumRows(mat);
  if (numRows != bfMatDenseComplexGetNumRows(otherMat))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize numCols = bfMatDenseComplexGetNumCols(mat);
  if (numCols != bfMatDenseComplexGetNumCols(otherMat))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize stride = op1->rowStride + op1->colStride;
  for (BfSize i = 0; i < op2->numElts; ++i)
    *(op1->data + i*stride) += op2->data[i];

  END_ERROR_HANDLING() {}
}

void bfMatDenseComplexCooComplexAddInplace(BfMatDenseComplex *op1,
                                           BfMatCooComplex const *op2) {
  BEGIN_ERROR_HANDLING();

  BfMat *mat = bfMatDenseComplexToMat(op1);
  BfMat const *otherMat = bfMatCooComplexConstToMatConst(op2);

  BfSize numRows = bfMatDenseComplexGetNumRows(mat);
  if (numRows != bfMatCooComplexGetNumRows(otherMat))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize numCols = bfMatDenseComplexGetNumCols(mat);
  if (numCols != bfMatCooComplexGetNumCols(otherMat))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  for (BfSize k = 0; k < op2->numElts; ++k) {
    BfSize i = op2->rowInd[k];
    BfSize j = op2->colInd[k];
    BfSize offset = i*op1->rowStride + j*op1->colStride;
    assert(offset < numRows*numCols);
    *(op1->data + offset) += op2->value[k];
  }

  END_ERROR_HANDLING() {}
}

void bfMatDenseComplexCooRealAddInplace(BfMatDenseComplex *op1,
                                        BfMatCooReal const *op2) {
  BEGIN_ERROR_HANDLING();

  BfMat *mat = bfMatDenseComplexToMat(op1);
  BfMat const *otherMat = bfMatCooRealConstToMatConst(op2);

  BfSize numRows = bfMatGetNumRows(mat);
  if (numRows != bfMatGetNumRows(otherMat))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize numCols = bfMatGetNumCols(mat);
  if (numCols != bfMatGetNumCols(otherMat))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  for (BfSize k = 0; k < op2->numElts; ++k) {
    BfSize i = op2->rowInd[k];
    BfSize j = op2->colInd[k];
    BfSize offset = i*op1->rowStride + j*op1->colStride;
    assert(offset < numRows*numCols);
    *(op1->data + offset) += op2->value[k];
  }

  END_ERROR_HANDLING() {}
}

BfMatDenseComplex *
bfMatDenseComplexDenseComplexSub(BfMatDenseComplex const *op1,
                                 BfMatDenseComplex const *op2)
{
  BEGIN_ERROR_HANDLING();

  BfMat const *mat1 = bfMatDenseComplexConstToMatConst(op1);
  if (bfMatIsTransposed(mat1))
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);

  BfMat const *mat2 = bfMatDenseComplexConstToMatConst(op2);
  if (bfMatIsTransposed(mat2))
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);

  BfSize m = bfMatGetNumRows(mat1);
  if (bfMatGetNumRows(mat2) != m)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize n = bfMatGetNumCols(mat2);
  if (bfMatGetNumCols(mat2) != n)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfMatDenseComplex *result = bfMatDenseComplexNew();
  HANDLE_ERROR();

  bfMatDenseComplexInit(result, m, n);
  HANDLE_ERROR();

  BfComplex *resultPtr = result->data;
  for (BfSize i = 0; i < m; ++i) {
    BfComplex const *rowPtr1 = op1->data + i*op1->rowStride;
    BfComplex const *rowPtr2 = op2->data + i*op2->rowStride;
    for (BfSize j = 0; j < n; ++j) {
      *resultPtr = *rowPtr1 - *rowPtr2;
      resultPtr += result->colStride;
    }
  }

  END_ERROR_HANDLING()
    bfMatDenseComplexDeinitAndDealloc(&result);

  return result;
}

BfMatDenseComplex *
bfMatDenseComplexDenseComplexMul(BfMatDenseComplex const *op1,
                                 BfMatDenseComplex const *op2)
{
  BEGIN_ERROR_HANDLING();

  enum CBLAS_TRANSPOSE transa = getCblasTranspose(op1);
  enum CBLAS_TRANSPOSE transb = getCblasTranspose(op2);

  BfMat const *super1 = bfMatDenseComplexConstToMatConst(op1);
  BfMat const *super2 = bfMatDenseComplexConstToMatConst(op2);

  BfMatDenseComplex *result = NULL;

  BfSize k = bfMatGetNumCols(super1);
  if (k != bfMatGetNumRows(super2))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  result = bfMatDenseComplexNew();
  HANDLE_ERROR();

  BfSize m = bfMatGetNumRows(super1);
  BfSize n = bfMatGetNumCols(super2);
  bfMatDenseComplexInit(result, m, n);
  HANDLE_ERROR();

  BfComplex alpha = 1, beta = 0;

  BfSize lda = getLeadingDimension(op1);
  BfSize ldb = getLeadingDimension(op2);
  BfSize ldc = getLeadingDimension(result);

  assert(m > 0 && n > 0 && k > 0);

  if (transa == CblasNoTrans) {
    assert(lda >= k);
  } else {
    assert(transa == CblasTrans || transa == CblasConjTrans);
    assert(lda >= m);
  }

  if (transb == CblasNoTrans) {
    assert(ldb >= n);
  } else {
    assert(transb == CblasTrans || transb == CblasConjTrans);
    assert(ldb >= k);
  }

  assert(ldc >= n);

  cblas_zgemm(CblasRowMajor, transa, transb, m, n, k,
              &alpha, op1->data, lda, op2->data, ldb, &beta, result->data, ldc);

  /* TODO: handle cblas errors  */

  END_ERROR_HANDLING() {
    if (result != NULL)
      bfMatDenseComplexDeinitAndDealloc(&result);
  }

  return result;
}

BfMatDenseComplex *
bfMatDenseComplexDenseComplexLstSq(BfMatDenseComplex const *lhs,
                                   BfMatDenseComplex const *rhs)
{
  BfMat const *lhsSuper = bfMatDenseComplexConstToMatConst(lhs);
  assert(!bfMatIsTransposed(lhsSuper));

  BfSize m = lhsSuper->numRows;
  BfSize n = lhsSuper->numCols;
  BfSize p = m < n ? m : n;

  BEGIN_ERROR_HANDLING();

  /* compute SVD of A */

  BfMatDenseComplex *U = bfMatDenseComplexNew();
  bfMatDenseComplexInit(U, m, p);

  BfMatDiagReal *S = bfMatDiagRealNew();
  bfMatDiagRealInit(S, p, p);

  BfMatDenseComplex *VH = bfMatDenseComplexNew();
  bfMatDenseComplexInit(VH, p, n);

  bfMatDenseComplexSvd(lhs, U, S, VH);
  HANDLE_ERROR();

  /* compute tolerance and compute number of terms to
   * retain in pseudoinverse */

  BfReal const atol = BF_EPS_MACH;
  BfReal const rtol = (m > n ? m : n)*BF_EPS_MACH;

  BfReal const *sigma = S->data;

  BfReal tol = rtol*sigma[0] + atol;

  BfSize k;
  for (k = 0; k < p; ++k)
    if (sigma[k] < tol)
      break;

  /* TODO: the casting below got kind of out of hand... clean it up */

  /* get subblocks of truncated SVD */

  BfMat *UkH = bfMatDenseComplexGetColRange(bfMatDenseComplexToMat(U), 0, k);
  bfMatConjTrans(UkH);

  BfMatDiagReal *Sk = bfMatDiagRealGetDiagBlock(S, 0, k);

  BfMat *Vk = bfMatGetRowRange(bfMatDenseComplexToMat(VH), 0, k);
  bfMatConjTrans(Vk);

  /* solve least squares problem */

  BfMat *tmp1 = bfMatMul(UkH, bfMatDenseComplexConstToMatConst(rhs));
  HANDLE_ERROR();

  BfMat *tmp2 = bfMatDenseComplexToMat(
    bfMatDiagRealDenseComplexSolve(Sk, bfMatToMatDenseComplex(tmp1)));
  HANDLE_ERROR();

  BfMat *result = bfMatMul(Vk, tmp2);
  HANDLE_ERROR();

  END_ERROR_HANDLING()
    bfMatDenseComplexDelete(&result);

  bfMatDenseComplexDeinitAndDealloc(&U);
  bfMatDiagRealDeinitAndDealloc(&S);
  bfMatDenseComplexDeinitAndDealloc(&VH);

  bfMatDenseComplexDelete(&UkH);
  bfMatDiagRealDeinitAndDealloc(&Sk);
  bfMatDenseComplexDelete(&Vk);

  bfMatDelete(&tmp1);
  bfMatDelete(&tmp2);

  return bfMatToMatDenseComplex(result);
}

bool bfMatDenseComplexIsFinite(BfMatDenseComplex const *mat) {
  for (BfSize i = 0; i < mat->super.numRows; ++i) {
    BfComplex *rowPtr = mat->data + i*mat->rowStride;
    for (BfSize j = 0; j < mat->super.numCols; ++j) {
      BfComplex elt = *(rowPtr + j*mat->colStride);
      if (!isfinite(creal(elt)) || !isfinite(cimag(elt)))
        return false;
    }
  }
  return true;
}

void bf_zmat_add_diag(BfMatDenseComplex *mat, BfReal value) {
  BfSize m = mat->super.numRows;
  BfSize n = mat->super.numCols;
  BfSize p = m < n ? m : n;
  for (BfSize i = 0; i < p; ++i)
    *(mat->data + i*mat->rowStride + i*mat->colStride) += value;
}

BfMatDenseComplex *
bfMatDenseComplexDenseComplexSolve(BfMatDenseComplex const *A,
                                   BfMatDenseComplex const *B) {
  BfSize m = A->super.numRows;
  if (B->super.numRows != m) {
    bfSetError(BF_ERROR_INVALID_ARGUMENTS);
    return NULL;
  }

  BfSize n = A->super.numCols;
  if (m != n) {
    bfSetError(BF_ERROR_INVALID_ARGUMENTS);
    return NULL;
  }

  BfSize p = B->super.numCols; /* == number of RHSs */

  BEGIN_ERROR_HANDLING();

  int *ipiv = malloc(m*sizeof(int));

  /* Create a copy of B here since X will be overwritten by LAPACK */
  BfMatDenseComplex *X = bfMatDenseComplexNew();
  bfMatDenseComplexInit(X, n, p);
  bfMatDenseComplexSet(X, B);

  lapack_int error = LAPACKE_zgesv(
    LAPACK_ROW_MAJOR, m, p, A->data, A->rowStride,
    ipiv, X->data, X->rowStride);

  if (error == LAPACK_WORK_MEMORY_ERROR)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  if (error != 0)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR)

  END_ERROR_HANDLING()
    free(ipiv);

  return X;
}

/** Upcasting: */

BfMat *bfMatDenseComplexToMat(BfMatDenseComplex *matDenseComplex) {
  return &matDenseComplex->super;
}

BfMat const *bfMatDenseComplexConstToMatConst(BfMatDenseComplex const *matDenseComplex) {
  return &matDenseComplex->super;
}

/** Downcasting: */

BfMatDenseComplex *bfMatToMatDenseComplex(BfMat *mat) {
  if (!bfMatInstanceOf(mat, BF_TYPE_MAT_DENSE_COMPLEX)) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return NULL;
  } else {
    return (BfMatDenseComplex *)mat;
  }
}

BfMatDenseComplex const *bfMatConstToMatDenseComplexConst(BfMat const *mat) {
  if (!bfMatInstanceOf(mat, BF_TYPE_MAT_DENSE_COMPLEX)) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return NULL;
  } else {
    return (BfMatDenseComplex const *)mat;
  }
}

/** Implementation: MatDenseComplex */

BfMatDenseComplex *bfMatDenseComplexNew() {
  BEGIN_ERROR_HANDLING();

  BfMatDenseComplex *mat = malloc(sizeof(BfMatDenseComplex));
  if (mat == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING() {}

  return mat;
}

BfMatDenseComplex *bfMatDenseComplexNewViewFromPtr(BfSize numRows, BfSize numCols, BfComplex *data) {
  BEGIN_ERROR_HANDLING();

  BfMatDenseComplex *matDenseComplex = bfMatDenseComplexNew();
  HANDLE_ERROR();

  bfMatDenseComplexInitViewFromPtr(matDenseComplex, numRows, numCols, data);

  END_ERROR_HANDLING() {
    assert(false);
  }

  return matDenseComplex;
}

BfMatDenseComplex *bfMatDenseComplexZeros(BfSize numRows, BfSize numCols) {
  BEGIN_ERROR_HANDLING();

  BfMatDenseComplex *zeros = NULL;

  zeros = bfMatDenseComplexNew();
  HANDLE_ERROR();

  bfMatDenseComplexInit(zeros, numRows, numCols);
  HANDLE_ERROR();

  for (BfSize i = 0; i < numRows*numCols; ++i)
    zeros->data[i] = 0;

  END_ERROR_HANDLING()
    bfMatDenseComplexDeinitAndDealloc(&zeros);

  return zeros;
}

BfMatDenseComplex *bfMatDenseComplexFromFile(char const *path, BfSize numRows, BfSize numCols) {
  BEGIN_ERROR_HANDLING();

  BfMatDenseComplex *matDenseComplex = bfMatDenseComplexNew();
  HANDLE_ERROR();

  bfMatDenseComplexInit(matDenseComplex, numRows, numCols);
  HANDLE_ERROR();

  FILE *fp = fopen(path, "r");
  if (fp == NULL)
    RAISE_ERROR(BF_ERROR_FILE_ERROR);

  BfSize size = numRows*numCols;
  if (fread(matDenseComplex->data, sizeof(BfComplex), size, fp) != size)
    RAISE_ERROR(BF_ERROR_FILE_ERROR);

  END_ERROR_HANDLING()
    bfMatDenseComplexDeinitAndDealloc(&matDenseComplex);

  return matDenseComplex;
}

void bfMatDenseComplexInit(BfMatDenseComplex *mat,
                           BfSize numRows, BfSize numCols) {
  BEGIN_ERROR_HANDLING();

  bfMatInit(&mat->super, &MAT_VTABLE, numRows, numCols);
  HANDLE_ERROR();

  mat->rowStride = numCols;
  mat->colStride = 1;

  mat->data = malloc(numRows*numCols*sizeof(BfComplex));
  if (mat->data == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING()
    bfMatDeinit(&mat->super);
}

void bfMatDenseComplexInitViewFromPtr(BfMatDenseComplex *matDenseComplex, BfSize numRows, BfSize numCols, BfComplex *data) {
  BEGIN_ERROR_HANDLING();

  bfMatInit(&matDenseComplex->super, &MAT_VTABLE, numRows, numCols);
  HANDLE_ERROR();

  matDenseComplex->super.props |= BF_MAT_PROPS_VIEW;

  matDenseComplex->rowStride = numCols;
  matDenseComplex->colStride = 1;
  matDenseComplex->data = data;

  END_ERROR_HANDLING() {
    assert(false);
  }
}

void bfMatDenseComplexDeinit(BfMatDenseComplex *mat) {
  if (!(mat->super.props & BF_MAT_PROPS_VIEW))
    free(mat->data);

  mat->data = NULL;
}

void bfMatDenseComplexDealloc(BfMatDenseComplex **mat) {
  free(*mat);
  *mat = NULL;
}

void bfMatDenseComplexDeinitAndDealloc(BfMatDenseComplex **mat) {
  bfMatDenseComplexDeinit(*mat);
  bfMatDenseComplexDealloc(mat);
}
