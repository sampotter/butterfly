#include <bf/mat_dense_complex.h>

#include <math.h>

#include <openblas/lapacke.h>

#include <bf/assert.h>
#include <bf/blas.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/lu_dense_complex.h>
#include <bf/macros.h>
#include <bf/mat_coo_complex.h>
#include <bf/mat_coo_real.h>
#include <bf/mat_dense_real.h>
#include <bf/mat_diag_real.h>
#include <bf/mem.h>
#include <bf/rand.h>
#include <bf/vec_complex.h>
#include <bf/vec_real.h>

#define NO_IMPORT_ARRAY
#include "numpy.h"

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
  BF_ERROR_BEGIN();

  BfVecReal *result = NULL;

  BfSize numRows = bfMatDenseComplexGetNumRows(matDenseComplex);
  if (numRows != bfMatDenseComplexGetNumRows(otherMatDenseComplex))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize numCols = bfMatDenseComplexGetNumCols(matDenseComplex);
  if (numCols != bfMatDenseComplexGetNumCols(otherMatDenseComplex))
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

  BF_ERROR_END()
    result = NULL;

  return bfVecRealToVec(result);
}

static BfVec *bfMatDenseComplexDenseComplexColDots(
  BfMatDenseComplex const *matDenseComplex,
  BfMatDenseComplex const *otherMatDenseComplex)
{
  BF_ERROR_BEGIN();

  BfVecComplex *result = NULL;

  BfSize numRows = bfMatDenseComplexGetNumRows(matDenseComplex);
  if (numRows != bfMatDenseComplexGetNumRows(otherMatDenseComplex))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize numCols = bfMatDenseComplexGetNumCols(matDenseComplex);
  if (numCols != bfMatDenseComplexGetNumCols(otherMatDenseComplex))
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

    BfComplex dot;
    cblas_zdotc_sub(
      /* n: */ numRows,
      /* x: */ colPtr,
      /* incx: */ matDenseComplex->rowStride,
      /* y: */ otherColPtr,
      /* incy: */ otherMatDenseComplex->rowStride,
      /* dotc: */ &dot);

    *resultPtr++ = dot;
  }

  BF_ERROR_END()
    result = NULL;

  return bfVecComplexToVec(result);
}

void bfMatDenseComplexScaleCols_real(BfMatDenseComplex *mat, BfVec const *vec) {
  BF_ERROR_BEGIN();

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

  BF_ERROR_END() {}
}

void bfMatDenseComplexScaleCols_complex(BfMatDenseComplex *mat, BfVec const *vec) {
  BF_ERROR_BEGIN();

  BfVecComplex const *vecComplex = bfVecConstToVecComplexConst(vec);
  HANDLE_ERROR();

  BfSize m = mat->super.numRows;
  BfSize n = mat->super.numCols;

  for (BfSize j = 0; j < n; ++j) {
    BfComplex const *readPtr = vecComplex->data + j*vecComplex->stride;
    BfComplex *writePtr = mat->data + j*mat->colStride;
    cblas_zscal(m, readPtr, writePtr, mat->rowStride);
  }

  BF_ERROR_END() {}
}

/** Interface: Mat */

static BfMatVtable MAT_VTABLE = {
  .GetView = (__typeof__(&bfMatDenseComplexGetView))bfMatDenseComplexGetView,
  .Copy = (__typeof__(&bfMatCopy))bfMatDenseComplexCopy,
  .Steal = (__typeof__(&bfMatSteal))bfMatDenseComplexSteal,
  .GetRowCopy = (__typeof__(&bfMatDenseComplexGetRowCopy))bfMatDenseComplexGetRowCopy,
  .GetRowView = (__typeof__(&bfMatDenseComplexGetRowView))bfMatDenseComplexGetRowView,
  .GetColView = (__typeof__(&bfMatDenseComplexGetColView))bfMatDenseComplexGetColView,
  .GetColRangeView = (__typeof__(&bfMatDenseComplexGetColRangeView))bfMatDenseComplexGetColRangeView,
  .Delete = (__typeof__(&bfMatDelete))bfMatDenseComplexDelete,
  .EmptyLike = (__typeof__(&bfMatDenseComplexEmptyLike))bfMatDenseComplexEmptyLike,
  .ZerosLike = (__typeof__(&bfMatDenseComplexZerosLike))bfMatDenseComplexZerosLike,
  .GetType = (__typeof__(&bfMatDenseComplexGetType))bfMatDenseComplexGetType,
  .NumBytes = (__typeof__(&bfMatNumBytes))bfMatDenseComplexNumBytes,
  .Save = (__typeof__(&bfMatDenseComplexSave))bfMatDenseComplexSave,
  .Print = (__typeof__(&bfMatDenseComplexPrint))bfMatDenseComplexPrint,
  .GetNumRows = (__typeof__(&bfMatGetNumRows))bfMatDenseComplexGetNumRows,
  .GetNumCols = (__typeof__(&bfMatGetNumCols))bfMatDenseComplexGetNumCols,
  .SetRow = (__typeof__(&bfMatDenseComplexSetRow))bfMatDenseComplexSetRow,
  .SetCol = (__typeof__(&bfMatDenseComplexSetCol))bfMatDenseComplexSetCol,
  .SetColRange = (__typeof__(&bfMatDenseComplexSetColRange))bfMatDenseComplexSetColRange,
  .GetRowRange = (__typeof__(&bfMatDenseComplexGetRowRange))bfMatDenseComplexGetRowRange,
  .GetRowRangeCopy = (__typeof__(&bfMatGetRowRangeCopy))bfMatDenseComplexGetRowRangeCopy,
  .GetColRange = (__typeof__(&bfMatGetColRange))bfMatDenseComplexGetColRange,
  .GetColRangeConst = (__typeof__(&bfMatGetColRangeConst))bfMatDenseComplexGetColRangeConst,
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
  .Rmul = (__typeof__(&bfMatRmul))bfMatDenseComplexRmul,
  .Solve = (__typeof__(&bfMatSolve))bfMatDenseComplexSolve,
  .SolveLU = (__typeof__(&bfMatDenseComplexSolveLU))bfMatDenseComplexSolveLU,
  .LstSq = (__typeof__(&bfMatDenseComplexLstSq))bfMatDenseComplexLstSq,
  .IsUpperTri = (__typeof__(&bfMatDenseComplexIsUpperTri))bfMatDenseComplexIsUpperTri,
  .BackwardSolveVec = (__typeof__(&bfMatDenseComplexBackwardSolveVec))bfMatDenseComplexBackwardSolveVec,
  .Negate = (__typeof__(&bfMatDenseComplexNegate))bfMatDenseComplexNegate,
  .ToType = (__typeof__(&bfMatToType))bfMatDenseComplexToType,
  .PrintBlocksDeep = (__typeof__(&bfMatPrintBlocksDeep))bfMatDenseComplexPrintBlocksDeep,
  .GetBlockView = (__typeof__(&bfMatGetBlockView))bfMatDenseComplexGetBlockView,
  .GetLu = (__typeof__(&bfMatGetLu))bfMatDenseComplexGetLu,
  .DivideCols = (__typeof__(&bfMatDivideCols))bfMatDenseComplexDivideCols,
  .Transpose = (__typeof__(&bfMatTranspose))bfMatDenseComplexTranspose,
};

BfMat *bfMatDenseComplexGetView(BfMat *mat) {
  BF_ERROR_BEGIN();

  BfMatDenseComplex *matDenseComplex = bfMatToMatDenseComplex(mat);
  HANDLE_ERROR();

  BfMatDenseComplex *matDenseComplexView = bfMatDenseComplexNew();

  *matDenseComplexView = *matDenseComplex;

  BfMat *matView = bfMatDenseComplexToMat(matDenseComplexView);

  matView->props |= BF_MAT_PROPS_VIEW;

  BF_ERROR_END()
    matView = NULL;

  return matView;
}

BfMat *bfMatDenseComplexCopy(BfMatDenseComplex const *matDenseComplex) {
  BF_ERROR_BEGIN();

  BfMatDenseComplex *copy = NULL;

  BfSize m = bfMatDenseComplexGetNumRows(matDenseComplex);
  BfSize n = bfMatDenseComplexGetNumCols(matDenseComplex);

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

  BF_ERROR_END()
    bfMatDenseComplexDeinitAndDealloc(&copy);

  return bfMatDenseComplexToMat(copy);
}

BfMat *bfMatDenseComplexSteal(BfMatDenseComplex *matDenseComplex) {
  BF_ERROR_BEGIN();

  BfMat *mat = bfMatDenseComplexToMat(matDenseComplex);

  if (bfMatIsView(mat))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfMatDenseComplex *matDenseComplexNew = bfMatDenseComplexNew();
  HANDLE_ERROR();

  *matDenseComplexNew = *matDenseComplex;

  mat->props |= BF_MAT_PROPS_VIEW;

  BF_ERROR_END() {
    BF_DIE();
  }

  return bfMatDenseComplexToMat(matDenseComplexNew);
}

BfVec *bfMatDenseComplexGetRowCopy(BfMat const *mat, BfSize i) {
  BF_ERROR_BEGIN();

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

  BF_ERROR_END() {}

  return bfVecComplexToVec(rowCopy);
}

BfVec *bfMatDenseComplexGetRowView(BfMat *mat, BfSize i) {
  BF_ERROR_BEGIN();

  BfMatDenseComplex *matDenseComplex = bfMatToMatDenseComplex(mat);
  HANDLE_ERROR();

  BfVecComplex *rowView = bfVecComplexNew();
  HANDLE_ERROR();

  BfSize n = bfMatGetNumCols(mat);
  BfSize stride = matDenseComplex->rowStride;
  BfComplex *data = matDenseComplex->data + i*matDenseComplex->rowStride;

  bfVecComplexInitView(rowView, n, stride, data);

  BF_ERROR_END()
    bfVecComplexDeinitAndDealloc(&rowView);

  return bfVecComplexToVec(rowView);
}

BfVec *bfMatDenseComplexGetColView(BfMat *mat, BfSize j) {
  BF_ERROR_BEGIN();

  BfMatDenseComplex *matDenseComplex = bfMatToMatDenseComplex(mat);
  HANDLE_ERROR();

  BfVecComplex *colView = bfVecComplexNew();
  HANDLE_ERROR();

  BfSize m = bfMatGetNumRows(mat);
  BfSize stride = matDenseComplex->rowStride;
  BfComplex *data = matDenseComplex->data + j*matDenseComplex->colStride;

  bfVecComplexInitView(colView, m, stride, data);

  BF_ERROR_END()
    bfVecComplexDeinitAndDealloc(&colView);

  return bfVecComplexToVec(colView);
}

BfVec *bfMatDenseComplexGetColRangeView(BfMat *mat, BfSize i0, BfSize i1, BfSize j) {
  BF_ERROR_BEGIN();

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

  BF_ERROR_END()
    bfVecComplexDeinitAndDealloc(&colView);

  return bfVecComplexToVec(colView);
}

void bfMatDenseComplexDelete(BfMatDenseComplex **matDenseComplex) {
  bfMatDenseComplexDeinitAndDealloc(matDenseComplex);
}

BfMat *bfMatDenseComplexEmptyLike(BfMat const *mat, BfSize numRows, BfSize numCols) {
  BF_ERROR_BEGIN();

  BfMatDenseComplex *result = NULL;

  if (numRows == BF_SIZE_BAD_VALUE)
    numRows = bfMatGetNumRows(mat);

  if (numCols == BF_SIZE_BAD_VALUE)
    numCols = bfMatGetNumCols(mat);

  result = bfMatDenseComplexNew();
  HANDLE_ERROR();

  bfMatDenseComplexInit(result, numRows, numCols);
  HANDLE_ERROR();

  BF_ERROR_END()
    bfMatDenseComplexDeinitAndDealloc(&result);

  return bfMatDenseComplexToMat(result);
}

BfMat *bfMatDenseComplexZerosLike(BfMat const *mat, BfSize numRows, BfSize numCols) {
  BF_ERROR_BEGIN();

  BfMatDenseComplex *zeros = NULL;

  if (numRows == BF_SIZE_BAD_VALUE)
    numRows = bfMatGetNumRows(mat);

  if (numCols == BF_SIZE_BAD_VALUE)
    numCols = bfMatGetNumCols(mat);

  zeros = bfMatDenseComplexZeros(numRows, numCols);
  HANDLE_ERROR();

  BF_ERROR_END()
    bfMatDenseComplexDeinitAndDealloc(&zeros);

  return bfMatDenseComplexToMat(zeros);
}

BfType bfMatDenseComplexGetType(BfMat const *mat) {
  (void)mat;
  return BF_TYPE_MAT_DENSE_COMPLEX;
}

BfSize bfMatDenseComplexNumBytes(BfMatDenseComplex const *matDenseComplex) {
  BfMat const *mat = bfMatDenseComplexConstToMatConst(matDenseComplex);
  return sizeof(BfComplex)*bfMatGetNumRows(mat)*bfMatGetNumCols(mat);
}

void bfMatDenseComplexSave(BfMat const *mat, char const *path) {
  BF_ERROR_BEGIN();

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

  BF_ERROR_END() {
    BF_DIE();
  }

  fclose(fp);
}

void bfMatDenseComplexPrint(BfMat const *mat, FILE *fp) {
  BF_ERROR_BEGIN();

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

  BF_ERROR_END() {}
}

BfSize bfMatDenseComplexGetNumRows(BfMatDenseComplex const *matDenseComplex) {
  BfMat const *mat = bfMatDenseComplexConstToMatConst(matDenseComplex);
  return bfMatIsTransposed(mat) ? mat->numCols : mat->numRows;
}

BfSize bfMatDenseComplexGetNumCols(BfMatDenseComplex const *matDenseComplex) {
  BfMat const *mat = bfMatDenseComplexConstToMatConst(matDenseComplex);
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
  BF_ERROR_BEGIN();

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

  BF_ERROR_END() {}
}

void bfMatDenseComplexSetCol(BfMat *mat, BfSize j, BfVec const *vec) {
  BF_ERROR_BEGIN();

  BfMatDenseComplex *matDenseComplex = NULL;
  BfVecComplex const *vecComplex = NULL;

  BfSize m = bfMatGetNumRows(mat);
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

  BF_ERROR_END() {}
}

void bfMatDenseComplexSetColRange(BfMat *mat, BfSize j, BfSize i0, BfSize i1,
                                  BfVec const *vec) {
  BF_ERROR_BEGIN();

  BfMatDenseComplex *matDenseComplex = NULL;
  BfVecComplex const *vecComplex = NULL;

  if (i0 > i1)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize m = bfMatGetNumRows(mat);
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

  BF_ERROR_END() {}
}

BfMat *bfMatDenseComplexGetRowRange(BfMat *mat, BfSize i0, BfSize i1) {
  BfSize numRows = mat->numRows;

  BF_ASSERT(i0 < i1);
  BF_ASSERT(i1 <= numRows);
  BF_ASSERT(!bfMatIsTransposed(mat)); // TODO: implement

  BF_ERROR_BEGIN();

  BfMat *matView = bfMatGetView(mat);

  BfMatDenseComplex *submat = bfMatToMatDenseComplex(matView);
  HANDLE_ERROR();

  if (i1 - i0 != numRows) {
    bfMatDenseComplexToMat(submat)->numRows = i1 - i0;
    submat->data += submat->rowStride*i0;
  }

  BF_ERROR_END()
    bfMatDelete(&matView);

  return matView;
}

BfMat *bfMatDenseComplexGetRowRangeCopy(BfMatDenseComplex const *matDenseComplex, BfSize i0, BfSize i1) {
  BF_ERROR_BEGIN();

  BfMat const *mat = bfMatDenseComplexConstToMatConst(matDenseComplex);
  HANDLE_ERROR();

  BfSize m = bfMatGetNumRows(mat);
  BfSize n = bfMatGetNumCols(mat);

  if (i0 >= i1)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  if (i1 > m)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  BfMatDenseComplex *matDenseComplexCopy = bfMatDenseComplexNew();
  HANDLE_ERROR();

  bfMatDenseComplexInit(matDenseComplexCopy, i1 - i0, n);
  HANDLE_ERROR();

  for (BfSize i = i0; i < i1; ++i) {
    BfComplex const *readPtr = matDenseComplex->data + i*matDenseComplex->rowStride;
    BfComplex *writePtr = matDenseComplexCopy->data + (i - i0)*matDenseComplexCopy->rowStride;
    for (BfSize j = 0; j < n; ++j) {
      *writePtr = *readPtr;
      readPtr += matDenseComplex->colStride;
      writePtr += matDenseComplexCopy->colStride;
    }
  }

  BF_ERROR_END() {
    bfMatDenseComplexDeinitAndDealloc(&matDenseComplexCopy);
  }

  return bfMatDenseComplexToMat(matDenseComplexCopy);
}

BfMat *bfMatDenseComplexGetColRange(BfMatDenseComplex *matDenseComplex, BfSize j0, BfSize j1) {
  BF_ERROR_BEGIN();

  BfMat *mat = bfMatDenseComplexToMat(matDenseComplex);

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

  BF_ERROR_END()
    bfMatDenseComplexDeinitAndDealloc(&matDenseComplexView);

  return bfMatDenseComplexToMat(matDenseComplexView);
}

BfMat const *bfMatDenseComplexGetColRangeConst(BfMatDenseComplex const *matDenseComplex, BfSize j0, BfSize j1) {
  return bfMatDenseComplexGetColRange((BfMatDenseComplex *)matDenseComplex, j0, j1);
}

void bfMatDenseComplexSetRowRange(BfMat *mat, BfSize i0, BfSize i1, BfMat const *rows) {
  BF_ERROR_BEGIN();

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

  BF_ERROR_END() {
    BF_DIE();
  }

  bfMatDenseComplexDelete(&matRows);
}

void bfMatDenseComplexPermuteRows(BfMat *mat, BfPerm const *perm) {
  BF_ERROR_BEGIN();

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

  BF_ERROR_END() {}

  bfMatDelete(&matPerm);
}

BfVec *bfMatDenseComplexColDists(BfMat const *mat, BfMat const *otherMat) {
  BF_ERROR_BEGIN();

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

  BF_ERROR_END()
    result = NULL;

  return result;
}

BfVec *bfMatDenseComplexColDots(BfMat const *mat, BfMat const *otherMat) {
  BF_ERROR_BEGIN();

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

  BF_ERROR_END()
    result = NULL;

  return result;
}

BfVec *bfMatDenseComplexColNorms(BfMat const *mat) {
  BF_ERROR_BEGIN();

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

  BF_ERROR_END()
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
  BF_ERROR_BEGIN();

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

  BF_ERROR_END() {}
}

void bfMatDenseComplexAddInplace(BfMat *mat, BfMat const *otherMat) {
  BF_ERROR_BEGIN();

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

  BF_ERROR_END() {}
}

BfMat *bfMatDenseComplexSub(BfMat const *mat, BfMat const *otherMat) {
  BF_ERROR_BEGIN();

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

  BF_ERROR_END() {}

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
  BF_ERROR_BEGIN();

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

  BF_ERROR_END() {}
}

static BfMat *mul_matDenseReal(BfMatDenseComplex const *op1, BfMatDenseReal const *op2) {
  BF_ERROR_BEGIN();

  BfSize m = bfMatDenseComplexGetNumRows(op1);
  BfSize n = bfMatDenseRealGetNumCols(op2);

  BfMatDenseComplex *res = bfMatDenseComplexNew();
  HANDLE_ERROR();

  bfMatDenseComplexInit(res, m, n);
  HANDLE_ERROR();

  /* BLAS doesn't specify mixed type GEMs, so we do it in two parts below.
   *
   * First, do "Re(res) = Re(op1)*op2": */

  BfMat const *op1Re = BF_TO_MAT(bfMatDenseComplexGetReViewConst(op1));
  HANDLE_ERROR();

  BfMat *resRe = BF_TO_MAT(bfMatDenseComplexGetReView(res));
  HANDLE_ERROR();

  bfMatMulInplace(op1Re, BF_TO_MAT(op2), resRe);
  HANDLE_ERROR();

  /* Now do "Im(res) = Im(op1)*op2": */

  BfMat const *op1Im = BF_TO_MAT(bfMatDenseComplexGetImViewConst(op1));
  HANDLE_ERROR();

  BfMat *resIm = BF_TO_MAT(bfMatDenseComplexGetImView(res));
  HANDLE_ERROR();

  bfMatMulInplace(op1Im, BF_TO_MAT(op2), resIm);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  return bfMatDenseComplexToMat(res);
}

BfMat *bfMatDenseComplexMul(BfMat const *op1, BfMat const *op2) {
  BF_ERROR_BEGIN();

  if (bfMatGetNumCols(op1) != bfMatGetNumRows(op2))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfMatDenseComplex const *matDenseComplex = bfMatConstToMatDenseComplexConst(op1);

  BfMat *result = NULL;

  switch (bfMatGetType(op2)) {
  case BF_TYPE_MAT_DENSE_COMPLEX:
    result = bfMatDenseComplexToMat(
      bfMatDenseComplexDenseComplexMul(
        matDenseComplex, bfMatConstToMatDenseComplexConst(op2)));
    break;
  case BF_TYPE_MAT_DENSE_REAL:
    result = mul_matDenseReal(matDenseComplex, bfMatConstToMatDenseRealConst(op2));
    break;
  default:
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);
  }

  BF_ERROR_END()
    bfMatDelete(&result);

  return result;
}

static BfVecComplex *
mulVec_complex(BfMatDenseComplex const *matDenseComplex,
               BfVecComplex const *vecComplex,
               BfSize m, BfSize n)
{
  BF_ERROR_BEGIN();

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

  BF_ERROR_END()
    bfVecComplexDeinitAndDealloc(&result);

  return result;
}

static BfMat *rmul_matDenseComplex(BfMatDenseComplex const *matDenseComplex, BfMatDenseComplex const *otherMatDenseComplex) {
  BF_ERROR_BEGIN();

  BfMatDenseComplex *result = NULL;

  BfSize m = bfMatDenseComplexGetNumRows(otherMatDenseComplex);
  BfSize n = bfMatDenseComplexGetNumCols(matDenseComplex);
  BfSize k = bfMatDenseComplexGetNumCols(otherMatDenseComplex);

  if (k != bfMatDenseComplexGetNumRows(matDenseComplex))
    RAISE_ERROR(BF_ERROR_INCOMPATIBLE_SHAPES);

  result = bfMatDenseComplexNew();
  HANDLE_ERROR();

  bfMatDenseComplexInit(result, m, n);
  HANDLE_ERROR();

  enum CBLAS_TRANSPOSE transa = getCblasTranspose(otherMatDenseComplex);
  enum CBLAS_TRANSPOSE transb = getCblasTranspose(matDenseComplex);

  BfComplex alpha = 1;
  BfComplex beta = 0;

  BfComplex *a = otherMatDenseComplex->data;
  BfComplex *b = matDenseComplex->data;
  BfComplex *c = result->data;

  BfSize lda = getLeadingDimension(otherMatDenseComplex);
  BfSize ldb = getLeadingDimension(matDenseComplex);
  BfSize ldc = getLeadingDimension(result);

  cblas_zgemm(CblasRowMajor, transa, transb, m, n, k, &alpha, a, lda, b, ldb, &beta, c, ldc);

  /* TODO: handle clbas errors... */

  BF_ERROR_END() {
    BF_DIE();
  }

  return bfMatDenseComplexToMat(result);
}

BfMat *bfMatDenseComplexRmul(BfMatDenseComplex const *matDenseComplex, BfMat const *mat) {
  switch (bfMatGetType(mat)) {
  case BF_TYPE_MAT_DENSE_COMPLEX:
    return rmul_matDenseComplex(matDenseComplex, bfMatConstToMatDenseComplexConst(mat));
  default:
    bfSetError(BF_ERROR_NOT_IMPLEMENTED);
    return NULL;
  }
}

BfVec *bfMatDenseComplexMulVec(BfMat const *mat, BfVec const *vec) {
  BF_ERROR_BEGIN();

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

  BF_ERROR_END()
    bfVecDelete(&result);

  return result;
}

static BfMat *
solve_matDenseComplex_tri(BfMatDenseComplex const *matDenseComplex, BfMat const *otherMat) {
  BF_ERROR_BEGIN();

  BfMat const *mat = bfMatDenseComplexConstToMatConst(matDenseComplex);

  BfSize m = bfMatGetNumRows(mat);
  BfSize n = bfMatGetNumCols(mat);
  BfSize p = bfMatGetNumCols(otherMat);

  if (m != n)
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);

  if (bfMatGetNumRows(otherMat) != m)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BF_ASSERT(mat->props & BF_MAT_PROPS_TRI);

  if (mat->props & BF_MAT_PROPS_TRANS)
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);

  BfMatDenseComplex const *otherMatDenseComplex = bfMatConstToMatDenseComplexConst(otherMat);

  BfMatDenseComplex *resultMatDenseComplex = bfMatDenseComplexNew();
  HANDLE_ERROR();

  bfMatDenseComplexInit(resultMatDenseComplex, n, p);
  HANDLE_ERROR();

  bfMemCopy(otherMatDenseComplex->data, n*p, sizeof(BfComplex), resultMatDenseComplex->data);

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

  BF_ERROR_END() {
    BF_DIE();
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
  BF_ERROR_BEGIN();

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

  BF_ERROR_END()
    bfMatDelete(&result);

  return result;
}

BfMat *bfMatDenseComplexLstSq(BfMat const *lhs, BfMat const *rhs) {
  BF_ERROR_BEGIN();

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

  BF_ERROR_END()
    bfMatDelete(&result);

  return result;
}

bool bfMatDenseComplexIsUpperTri(BfMat const *mat) {
  BF_ERROR_BEGIN();

  bool upperTri = true;

  BfMatDenseComplex const *matDenseComplex = NULL;

  matDenseComplex = bfMatConstToMatDenseComplexConst(mat);
  HANDLE_ERROR();

  BfSize m = bfMatGetNumRows(mat);

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

  BF_ERROR_END() {}

  return upperTri;
}

static BfVec *
backwardSolveVec_complex(BfMatDenseComplex const *matDenseComplex,
                         BfVecComplex const *vecComplex, BfSize m)
{
  BF_ERROR_BEGIN();

  BfVecComplex *result = bfVecToVecComplex(bfVecCopy(&vecComplex->super));
  HANDLE_ERROR();

  enum CBLAS_TRANSPOSE trans = getCblasTranspose(matDenseComplex);
  BfSize lda = getLeadingDimension(matDenseComplex);

  cblas_ztrsv(CblasRowMajor, CblasUpper, trans, CblasNonUnit, m,
              matDenseComplex->data, lda, result->data, result->stride);

  BF_ERROR_END() {}

  return bfVecComplexToVec(result);
}

BfVec *bfMatDenseComplexBackwardSolveVec(BfMat const *mat, BfVec const *vec) {
  BF_ERROR_BEGIN();

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

  BF_ERROR_END()
    bfVecDelete(&result);

  return result;
}

void bfMatDenseComplexNegate(BfMat *mat) {
  BF_ERROR_BEGIN();

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

  BF_ERROR_END() {}
}

BfMat *bfMatDenseComplexToType(BfMatDenseComplex const *matDenseComplex, BfType type) {
  switch (type) {
  case BF_TYPE_MAT_DENSE_COMPLEX:
    return bfMatDenseComplexCopy(matDenseComplex);
  default:
    bfSetError(BF_ERROR_NOT_IMPLEMENTED);
    return NULL;
  }
}

void bfMatDenseComplexPrintBlocksDeep(BfMatDenseComplex const *matDenseComplex, FILE *fp, BfSize i0, BfSize j0, BfSize depth) {
  BfMat const *mat = bfMatDenseComplexConstToMatConst(matDenseComplex);

  BfSize i1 = i0 + bfMatGetNumRows(mat);
  BfSize j1 = j0 + bfMatGetNumCols(mat);

  fprintf(fp, "%u %lu %lu %lu %lu %lu\n", BF_TYPE_MAT_DENSE_COMPLEX, i0, i1, j0, j1, depth);
}

BfMat *bfMatDenseComplexGetBlockView(BfMatDenseComplex *mat, BfSize i0, BfSize i1, BfSize j0, BfSize j1) {
  BF_ERROR_BEGIN();

  BfMatDenseComplex *blockView = bfMatDenseComplexNew();
  HANDLE_ERROR();

  bfMatInit(&blockView->super, &MAT_VTABLE, i1 - i0, j1 - j0);

  blockView->super.props |= BF_MAT_PROPS_VIEW;

  blockView->rowStride = mat->rowStride;
  blockView->colStride = mat->colStride;
  blockView->data = mat->data + i0*mat->rowStride + j0*mat->colStride;

  BF_ERROR_END() {}

  return bfMatDenseComplexToMat(blockView);
}

BfLu *bfMatDenseComplexGetLu(BfMatDenseComplex const *matDenseComplex) {
  BF_ERROR_BEGIN();

  BfMat const *mat = bfMatDenseComplexConstToMatConst(matDenseComplex);

  BfLuDenseComplex *luDenseComplex = bfLuDenseComplexNew();
  HANDLE_ERROR();

  bfLuDenseComplexInit(luDenseComplex, mat);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  return bfLuDenseComplexToLu(luDenseComplex);
}

static void divideCols_vecReal(BfMatDenseComplex *matDenseComplex, BfVec const *vec) {
  BF_ERROR_BEGIN();

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

  BF_ERROR_END() {
    BF_DIE();
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

void bfMatDenseComplexTranspose(BfMatDenseComplex *matDenseComplex) {
  BfMat *mat = bfMatDenseComplexToMat(matDenseComplex);
  bfMatConjTrans(mat);
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
  BF_ERROR_BEGIN();

  BfSize dstRows = bfMatDenseComplexGetNumRows(dst);
  if (dstRows != bfMatDenseComplexGetNumRows(src))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize dstCols = bfMatDenseComplexGetNumCols(dst);
  if (dstCols != bfMatDenseComplexGetNumCols(src))
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

  BF_ERROR_END() {}
}

void bfMatDenseComplexSvd(BfMatDenseComplex const *mat, BfMatDenseComplex *U,
                          BfMatDiagReal *S, BfMatDenseComplex *VH) {
  BF_ERROR_BEGIN();

  BfMat const *super = bfMatDenseComplexConstToMatConst(mat);

  BfSize m = super->numRows;
  BfSize n = super->numCols;

  BfReal *superb = NULL;

  /* zgesvd will overwrite A, so allocate space for a copy */
  BfComplex *dataCopy = bfMemAlloc(m*n, sizeof(BfComplex));
  if (dataCopy == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  /* copy contents of A */
  bfMemCopy(mat->data, m*n, sizeof(BfComplex), dataCopy);
  /* TODO: error-handling */

  /* output array which contains information about superdiagonal
   * elements which didn't converge
   *
   * more info here: tinyurl.com/2p8f5ev3 */
  superb = bfMemAlloc(((m < n) ? m : n) - 1, sizeof(BfReal));
  if (superb == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  /* compute the SVD */
  lapack_int info = LAPACKE_zgesvd(
    LAPACK_ROW_MAJOR, 'S', 'S', m, n, dataCopy, n, S->data,
    U->data, m < n ? m : n, VH->data, n, superb);

  /* check for invalid arguments */
  if (info < 0)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  /* check for errors */
  if (info > 0)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  BF_ERROR_END() {}

  bfMatDenseComplexToMat(U)->props |= BF_MAT_PROPS_ORTHO;
  bfMatDenseComplexToMat(VH)->props |= BF_MAT_PROPS_ORTHO;

  free(dataCopy);
  free(superb);
}

void bfMatDenseComplexDenseComplexAddInplace(BfMatDenseComplex *op1,
                                             BfMatDenseComplex const *op2)
{
  BF_ERROR_BEGIN();

  BfSize numRows = bfMatDenseComplexGetNumRows(op1);
  if (numRows != bfMatDenseComplexGetNumRows(op2))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize numCols = bfMatDenseComplexGetNumCols(op1);
  if (numCols != bfMatDenseComplexGetNumCols(op2))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  for (BfSize i = 0; i < numRows*numCols; ++i)
    op1->data[i] += op2->data[i];

  BF_ERROR_END() {}
}

void bfMatDenseComplexDiagRealAddInplace(BfMatDenseComplex *op1,
                                         BfMatDiagReal const *op2) {
  BF_ERROR_BEGIN();

  BfSize numRows = bfMatDenseComplexGetNumRows(op1);
  if (numRows != bfMatDiagRealGetNumRows(op2))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize numCols = bfMatDenseComplexGetNumCols(op1);
  if (numCols != bfMatDiagRealGetNumCols(op2))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize stride = op1->rowStride + op1->colStride;
  for (BfSize i = 0; i < op2->numElts; ++i)
    *(op1->data + i*stride) += op2->data[i];

  BF_ERROR_END() {}
}

void bfMatDenseComplexCooComplexAddInplace(BfMatDenseComplex *op1,
                                           BfMatCooComplex const *op2) {
  BF_ERROR_BEGIN();

  BfSize numRows = bfMatDenseComplexGetNumRows(op1);
  if (numRows != bfMatCooComplexGetNumRows(op2))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize numCols = bfMatDenseComplexGetNumCols(op1);
  if (numCols != bfMatCooComplexGetNumCols(op2))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  for (BfSize k = 0; k < op2->numElts; ++k) {
    BfSize i = op2->rowInd[k];
    BfSize j = op2->colInd[k];
    BfSize offset = i*op1->rowStride + j*op1->colStride;
    BF_ASSERT(offset < numRows*numCols);
    *(op1->data + offset) += op2->value[k];
  }

  BF_ERROR_END() {}
}

void bfMatDenseComplexCooRealAddInplace(BfMatDenseComplex *op1,
                                        BfMatCooReal const *op2) {
  BF_ERROR_BEGIN();

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
    BF_ASSERT(offset < numRows*numCols);
    *(op1->data + offset) += op2->value[k];
  }

  BF_ERROR_END() {}
}

BfMatDenseComplex *
bfMatDenseComplexDenseComplexSub(BfMatDenseComplex const *op1,
                                 BfMatDenseComplex const *op2)
{
  BF_ERROR_BEGIN();

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
      rowPtr1 += op1->colStride;
      rowPtr2 += op2->colStride;
    }
  }

  BF_ERROR_END()
    bfMatDenseComplexDeinitAndDealloc(&result);

  return result;
}

BfMatDenseComplex *
bfMatDenseComplexDenseComplexMul(BfMatDenseComplex const *op1,
                                 BfMatDenseComplex const *op2)
{
  BF_ERROR_BEGIN();

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

  BF_ASSERT(m > 0 && n > 0 && k > 0);

  if (transa == CblasNoTrans) {
    BF_ASSERT(lda >= k);
  } else {
    BF_ASSERT(transa == CblasTrans || transa == CblasConjTrans);
    BF_ASSERT(lda >= m);
  }

  if (transb == CblasNoTrans) {
    BF_ASSERT(ldb >= n);
  } else {
    BF_ASSERT(transb == CblasTrans || transb == CblasConjTrans);
    BF_ASSERT(ldb >= k);
  }

  BF_ASSERT(ldc >= n);

  cblas_zgemm(CblasRowMajor, transa, transb, m, n, k,
              &alpha, op1->data, lda, op2->data, ldb, &beta, result->data, ldc);

  /* TODO: handle cblas errors  */

  BF_ERROR_END() {
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
  BF_ASSERT(!bfMatIsTransposed(lhsSuper));

  BfSize m = lhsSuper->numRows;
  BfSize n = lhsSuper->numCols;
  BfSize p = m < n ? m : n;

  BF_ERROR_BEGIN();

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

  BfMat *UkH = bfMatDenseComplexGetColRange(U, 0, k);
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

  BF_ERROR_END() {
    BF_DIE();
  }

  bfMatDenseComplexDeinitAndDealloc(&U);
  bfMatDiagRealDeinitAndDealloc(&S);
  bfMatDenseComplexDeinitAndDealloc(&VH);

  bfMatDelete(&UkH);
  bfMatDiagRealDeinitAndDealloc(&Sk);
  bfMatDelete(&Vk);

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

  BF_ERROR_BEGIN();

  int *ipiv = bfMemAlloc(m, sizeof(int));

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

  BF_ERROR_END() {
    BF_DIE();
  }

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
  BF_ERROR_BEGIN();

  BfMatDenseComplex *mat = bfMemAlloc(1, sizeof(BfMatDenseComplex));
  if (mat == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BF_ERROR_END() {}

  return mat;
}

BfMatDenseComplex *bfMatDenseComplexNewViewFromPtr(BfSize numRows, BfSize numCols, BfComplex *data) {
  BF_ERROR_BEGIN();

  BfMatDenseComplex *matDenseComplex = bfMatDenseComplexNew();
  HANDLE_ERROR();

  bfMatDenseComplexInitViewFromPtr(matDenseComplex, numRows, numCols, data);

  BF_ERROR_END() {
    BF_DIE();
  }

  return matDenseComplex;
}

BfMatDenseComplex *bfMatDenseComplexNewRandn(BfSize numRows, BfSize numCols) {
  BF_ERROR_BEGIN();

  BfMatDenseComplex *matDenseComplex = bfMemAlloc(1, sizeof(BfMatDenseComplex));
  HANDLE_ERROR();

  bfMatDenseComplexInitRandn(matDenseComplex, numRows, numCols);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  return matDenseComplex;
}

BfMatDenseComplex *bfMatDenseComplexNewViewFromPyArray(BfPtr *pyArray) {
  BF_ERROR_BEGIN();

  BfMatDenseComplex *matDenseComplex = bfMemAlloc(1, sizeof(BfMatDenseComplex));
  HANDLE_ERROR();

  bfMatDenseComplexInitViewFromPyArray(matDenseComplex, pyArray);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  return matDenseComplex;
}

BfMatDenseComplex *bfMatDenseComplexZeros(BfSize numRows, BfSize numCols) {
  BF_ERROR_BEGIN();

  BfMatDenseComplex *zeros = NULL;

  zeros = bfMatDenseComplexNew();
  HANDLE_ERROR();

  bfMatDenseComplexInit(zeros, numRows, numCols);
  HANDLE_ERROR();

  for (BfSize i = 0; i < numRows*numCols; ++i)
    zeros->data[i] = 0;

  BF_ERROR_END() {
    BF_DIE();
  }

  return zeros;
}

BfMatDenseComplex *bfMatDenseComplexFromFile(char const *path, BfSize numRows, BfSize numCols) {
  BF_ERROR_BEGIN();

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

  BF_ERROR_END()
    bfMatDenseComplexDeinitAndDealloc(&matDenseComplex);

  return matDenseComplex;
}

BfMatDenseComplex *bfMatDenseComplexNewIdentity(BfSize numRows, BfSize numCols) {
  BF_ERROR_BEGIN();

  BfMatDenseComplex *matDenseComplex = bfMatDenseComplexNew();
  HANDLE_ERROR();

  bfMatDenseComplexInitIdentity(matDenseComplex, numRows, numCols);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  return matDenseComplex;
}

void bfMatDenseComplexInit(BfMatDenseComplex *mat,
                           BfSize numRows, BfSize numCols) {
  BF_ERROR_BEGIN();

  bfMatInit(&mat->super, &MAT_VTABLE, numRows, numCols);
  HANDLE_ERROR();

  mat->rowStride = numCols;
  mat->colStride = 1;

  mat->data = bfMemAlloc(numRows*numCols, sizeof(BfComplex));
  if (mat->data == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  mat->pyArray = NULL;

  BF_ERROR_END()
    bfMatDeinit(&mat->super);
}

void bfMatDenseComplexInitViewFromPtr(BfMatDenseComplex *matDenseComplex, BfSize numRows, BfSize numCols, BfComplex *data) {
  BF_ERROR_BEGIN();

  bfMatInit(&matDenseComplex->super, &MAT_VTABLE, numRows, numCols);
  HANDLE_ERROR();

  matDenseComplex->super.props |= BF_MAT_PROPS_VIEW;

  matDenseComplex->rowStride = numCols;
  matDenseComplex->colStride = 1;
  matDenseComplex->data = data;

  matDenseComplex->pyArray = NULL;

  BF_ERROR_END() {
    BF_DIE();
  }
}

void bfMatDenseComplexInitRandn(BfMatDenseComplex *matDenseComplex, BfSize numRows, BfSize numCols) {
  BF_ERROR_BEGIN();

  bfMatDenseComplexInit(matDenseComplex, numRows, numCols);
  HANDLE_ERROR();

  bfComplexRandn(numRows*numCols, matDenseComplex->data);

  BF_ERROR_END() {
    BF_DIE();
  }
}

void bfMatDenseComplexInitViewFromPyArray(BfMatDenseComplex *matDenseComplex, BfPtr *pyArray) {
  BF_ERROR_BEGIN();

  PyObject *obj = (PyObject *)pyArray;
  if (!PyArray_Check(obj))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  PyArrayObject *arrayObj = (PyArrayObject *)obj;
  BfSize numRows = PyArray_DIM(arrayObj, 0);
  BfSize numCols = PyArray_DIM(arrayObj, 1);
  BfComplex *data = PyArray_DATA(arrayObj);

  bfMatDenseComplexInitViewFromPtr(matDenseComplex, numRows, numCols, data);
  HANDLE_ERROR();

  matDenseComplex->pyArray = pyArray;

  Py_INCREF(obj);

  BF_ERROR_END() {
    BF_DIE();
  }
}

void bfMatDenseComplexInitIdentity(BfMatDenseComplex *matDenseComplex, BfSize numRows, BfSize numCols) {
  BF_ERROR_BEGIN();

  bfMatDenseComplexInit(matDenseComplex, numRows, numCols);
  HANDLE_ERROR();

  for (BfSize i = 0; i < numRows*numCols; ++i)
    matDenseComplex->data[i] = 0;

  for (BfSize i = 0; i < numRows; ++i)
    matDenseComplex->data[numCols*i + i] = 1;

  BF_ERROR_END() {
    BF_DIE();
  }
}

void bfMatDenseComplexDeinit(BfMatDenseComplex *mat) {
  if (!(mat->super.props & BF_MAT_PROPS_VIEW))
    free(mat->data);

  mat->data = NULL;

  // if (mat->pyArray != NULL) {
  //   PyObject *obj = (PyObject *)mat->pyArray;
  //   // Py_DECREF(obj);
  //   mat->pyArray = NULL;
  // }

  bfMatDeinit(&mat->super);
}

void bfMatDenseComplexDealloc(BfMatDenseComplex **mat) {
  free(*mat);
  *mat = NULL;
}

void bfMatDenseComplexDeinitAndDealloc(BfMatDenseComplex **mat) {
  bfMatDenseComplexDeinit(*mat);
  bfMatDenseComplexDealloc(mat);
}

BfComplex *bfMatDenseComplexGetDataPtr(BfMatDenseComplex *matDenseComplex) {
  return matDenseComplex->data;
}

BfComplex const *bfMatDenseComplexGetDataConstPtr(BfMatDenseComplex const *matDenseComplex) {
  return matDenseComplex->data;
}

static void setBlock_denseComplex(BfMatDenseComplex *matDenseComplex, BfSize i0, BfSize j0, BfMatDenseComplex const *blockMatDenseComplex) {
  BF_ERROR_BEGIN();

  BfMat *mat = bfMatDenseComplexToMat(matDenseComplex);
  BfSize numRows = bfMatGetNumRows(mat);
  BfSize numCols = bfMatGetNumCols(mat);

  BfMat const *blockMat = bfMatDenseComplexConstToMatConst(blockMatDenseComplex);
  BfSize numBlockRows = bfMatGetNumRows(blockMat);
  BfSize numBlockCols = bfMatGetNumCols(blockMat);

  BfSize i1 = i0 + numBlockRows;
  BfSize j1 = j0 + numBlockCols;

  if (i1 > numRows)
    RAISE_ERROR(BF_ERROR_OUT_OF_RANGE);

  if (j1 > numCols)
    RAISE_ERROR(BF_ERROR_OUT_OF_RANGE);

  for (BfSize i = i0; i < i1; ++i) {
    BfComplex *writePtr = matDenseComplex->data +
      i*matDenseComplex->rowStride + j0*matDenseComplex->colStride;
    BfComplex const *readPtr =
      blockMatDenseComplex->data + (i - i0)*blockMatDenseComplex->rowStride;
    for (BfSize j = j0; j < j1; ++j) {
      *writePtr = *readPtr;
      writePtr += matDenseComplex->colStride;
      readPtr += blockMatDenseComplex->colStride;
    }
  }

  BF_ERROR_END() {
    BF_DIE();
  }
}

void bfMatDenseComplexSetBlock(BfMatDenseComplex *matDenseComplex, BfSize i0, BfSize j0, BfMat const *mat) {
  switch (bfMatGetType(mat)) {
  case BF_TYPE_MAT_DENSE_COMPLEX:
    setBlock_denseComplex(matDenseComplex, i0, j0, bfMatConstToMatDenseComplexConst(mat));
    break;
  default:
    bfSetError(BF_ERROR_NOT_IMPLEMENTED);
  }
}

BfMatDenseReal *bfMatDenseComplexGetReView(BfMatDenseComplex *matDenseComplex) {
  BF_ERROR_BEGIN();

  if (matDenseComplex->rowStride != bfMatDenseComplexGetNumCols(matDenseComplex))
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);

  if (matDenseComplex->colStride != 1)
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);

  BfSize numRows = bfMatDenseComplexGetNumRows(matDenseComplex);
  BfSize numCols = bfMatDenseComplexGetNumCols(matDenseComplex);
  BfReal *ptr = (BfReal *)matDenseComplex->data;

  BF_ASSERT(matDenseComplex->rowStride == bfMatDenseComplexGetNumCols(matDenseComplex));
  BF_ASSERT(matDenseComplex->colStride == 1);
  BfSize rowStride = 2*numCols;
  BfSize colStride = 2;

  BfMatDenseReal *reView = bfMatDenseRealNewViewFromPtr(numRows, numCols, ptr, rowStride, colStride);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  return reView;
}

BfMatDenseReal const *bfMatDenseComplexGetReViewConst(BfMatDenseComplex const *matDenseComplex) {
  BF_ERROR_BEGIN();

  if (matDenseComplex->rowStride != bfMatDenseComplexGetNumCols(matDenseComplex))
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);

  if (matDenseComplex->colStride != 1)
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);

  BfSize numRows = bfMatDenseComplexGetNumRows(matDenseComplex);
  BfSize numCols = bfMatDenseComplexGetNumCols(matDenseComplex);
  BfReal const *ptr = (BfReal const *)matDenseComplex->data;

  BF_ASSERT(matDenseComplex->rowStride == bfMatDenseComplexGetNumCols(matDenseComplex));
  BF_ASSERT(matDenseComplex->colStride == 1);
  BfSize rowStride = 2*numCols;
  BfSize colStride = 2;

  BfMatDenseReal const *reViewConst = bfMatDenseRealNewViewConstFromPtr(numRows, numCols, ptr, rowStride, colStride);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  return reViewConst;
}

BfMatDenseReal *bfMatDenseComplexGetImView(BfMatDenseComplex *matDenseComplex) {
  BF_ERROR_BEGIN();

  if (matDenseComplex->rowStride != bfMatDenseComplexGetNumCols(matDenseComplex))
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);

  if (matDenseComplex->colStride != 1)
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);

  BfSize numRows = bfMatDenseComplexGetNumRows(matDenseComplex);
  BfSize numCols = bfMatDenseComplexGetNumCols(matDenseComplex);
  BfReal *ptr = ((BfReal *)matDenseComplex->data) + 1;

  BF_ASSERT(matDenseComplex->rowStride == bfMatDenseComplexGetNumCols(matDenseComplex));
  BF_ASSERT(matDenseComplex->colStride == 1);
  BfSize rowStride = 2*numCols;
  BfSize colStride = 2;

  BfMatDenseReal *imView = bfMatDenseRealNewViewFromPtr(numRows, numCols, ptr, rowStride, colStride);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  return imView;
}

BfMatDenseReal const *bfMatDenseComplexGetImViewConst(BfMatDenseComplex const *matDenseComplex) {
  BF_ERROR_BEGIN();

  if (matDenseComplex->rowStride != bfMatDenseComplexGetNumCols(matDenseComplex))
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);

  if (matDenseComplex->colStride != 1)
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);

  BfSize numRows = bfMatDenseComplexGetNumRows(matDenseComplex);
  BfSize numCols = bfMatDenseComplexGetNumCols(matDenseComplex);
  BfReal const *ptr = ((BfReal const *)matDenseComplex->data) + 1;

  BF_ASSERT(matDenseComplex->rowStride == bfMatDenseComplexGetNumCols(matDenseComplex));
  BF_ASSERT(matDenseComplex->colStride == 1);
  BfSize rowStride = 2*numCols;
  BfSize colStride = 2;

  BfMatDenseReal const *imViewConst = bfMatDenseRealNewViewConstFromPtr(numRows, numCols, ptr, rowStride, colStride);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  return imViewConst;
}

BfComplex *bfMatDenseComplexGetRowPtr(BfMatDenseComplex *matDenseComplex, BfSize i) {
  return matDenseComplex->data + i*matDenseComplex->rowStride;
}

BfComplex const *bfMatDenseComplexGetRowConstPtr(BfMatDenseComplex const *matDenseComplex, BfSize i) {
  return matDenseComplex->data + i*matDenseComplex->rowStride;
}
