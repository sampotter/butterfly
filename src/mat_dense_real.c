#include <bf/mat_dense_real.h>

#include <string.h>

#include <bf/real_array.h>
#include <bf/assert.h>
#include <bf/blas.h>
#include <bf/const.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mat_block.h>
#include <bf/mat_diag_real.h>
#include <bf/mem.h>
#include <bf/vec_real.h>

#include <string.h>

#include "macros.h"

static enum CBLAS_TRANSPOSE getCblasTranspose(BfMatDenseReal const *mat) {
  BfMat const *super = bfMatDenseRealConstToMatConst(mat);
  if (super->props & (BF_MAT_PROPS_TRANS | BF_MAT_PROPS_CONJ))
    return CblasConjTrans;
  else if (super->props & BF_MAT_PROPS_TRANS)
    return CblasTrans;
  else
    return CblasNoTrans;
}

static BfSize getLeadingDimension(BfMatDenseReal const *mat) {
  return mat->super.rowStride;
}

/** Interface: Mat */

static BfMatVtable MAT_VTABLE = {
  .GetView = (__typeof__(&bfMatGetView))bfMatDenseRealGetView,
  .Copy = (__typeof__(&bfMatCopy))bfMatDenseRealCopy,
  .Steal = (__typeof__(&bfMatSteal))bfMatDenseRealSteal,
  .GetRowView = (__typeof__(&bfMatGetRowView))bfMatDenseRealGetRowView,
  .GetColView = (__typeof__(&bfMatGetColView))bfMatDenseRealGetColView,
  .Delete = (__typeof__(&bfMatDelete))bfMatDenseRealDelete,
  .ZerosLike = (__typeof__(&bfMatZerosLike))bfMatDenseRealZerosLike,
  .GetType = (__typeof__(&bfMatGetType))bfMatDenseRealGetType,
  .NumBytes = (__typeof__(&bfMatNumBytes))bfMatDenseRealNumBytes,
  .Save = (__typeof__(&bfMatSave))bfMatDenseRealSave,
  .Dump = (__typeof__(&bfMatDump))bfMatDenseRealDump,
  .Print = (__typeof__(&bfMatPrint))bfMatDenseRealPrint,
  .GetNumRows = (__typeof__(&bfMatGetNumRows))bfMatDenseRealGetNumRows,
  .GetNumCols = (__typeof__(&bfMatGetNumCols))bfMatDenseRealGetNumCols,
  .SetRow = (__typeof__(&bfMatSetRow))bfMatDenseRealSetRow,
  .SetCol = (__typeof__(&bfMatSetCol))bfMatDenseRealSetCol,
  .GetRowRange = (__typeof__(&bfMatGetRowRange))bfMatDenseRealGetRowRange,
  .GetRowRangeCopy = (__typeof__(&bfMatGetRowRangeCopy))bfMatDenseRealGetRowRangeCopy,
  .GetColRangeCopy = (__typeof__(&bfMatGetColRangeCopy))bfMatDenseRealGetColRangeCopy,
  .PermuteRows = (__typeof__(&bfMatPermuteRows))bfMatDenseRealPermuteRows,
  .ScaleRows = (__typeof__(&bfMatScaleRows))bfMatDenseRealScaleRows,
  .Mul = (__typeof__(&bfMatMul))bfMatDenseRealMul,
  .MulVec = (__typeof__(&bfMatMulVec))bfMatDenseRealMulVec,
  .MulInplace = (__typeof__(&bfMatMulInplace))bfMatDenseRealMulInplace,
  .RmulVec = (__typeof__(&bfMatRmulVec))bfMatDenseRealRmulVec,
  .PrintBlocksDeep = (__typeof__(&bfMatPrintBlocksDeep))bfMatDenseRealPrintBlocksDeep,
  .NormMax = (__typeof__(&bfMatNormMax))bfMatDenseRealNormMax,
  .DistMax = (__typeof__(&bfMatDistMax))bfMatDenseRealDistMax,
};

BfMat *bfMatDenseRealGetView(BfMat *mat) {
  BF_ERROR_BEGIN();

  BfMatDenseReal *matDenseReal = bfMatToMatDenseReal(mat);
  HANDLE_ERROR();

  BfMatDenseReal *matDenseRealView = bfMatDenseRealNew();

  *matDenseRealView = *matDenseReal;

  BfMat *matView = bfMatDenseRealToMat(matDenseRealView);

  matView->props |= BF_MAT_PROPS_VIEW;

  BF_ERROR_END()
    matView = NULL;

  return matView;
}

BfMat *bfMatDenseRealCopy(BfMatDenseReal const *matDenseReal) {
  BF_ERROR_BEGIN();

  BfMatDenseReal *matDenseRealCopy = bfMatDenseRealNew();
  HANDLE_ERROR();

  bfMatDenseRealInitCopy(matDenseRealCopy, matDenseReal);
  HANDLE_ERROR();

  BF_ERROR_END()
    bfMatDenseRealDeinitAndDealloc(&matDenseRealCopy);

  return bfMatDenseRealToMat(matDenseRealCopy);
}

BfMat *bfMatDenseRealSteal(BfMatDenseReal *matDenseReal) {
  BF_ERROR_BEGIN();

  BfMat *mat = bfMatDenseRealToMat(matDenseReal);

  if (bfMatIsView(mat))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfMatDenseReal *matDenseRealNew = bfMatDenseRealNew();
  HANDLE_ERROR();

  *matDenseRealNew = *matDenseReal;

  mat->props |= BF_MAT_PROPS_VIEW;

  BF_ERROR_END() {
    BF_DIE();
  }

  return bfMatDenseRealToMat(matDenseRealNew);
}

BfVecReal *bfMatDenseRealGetRowView(BfMatDenseReal *matDenseReal, BfSize i) {
  BF_ERROR_BEGIN();

  BfMatDense *matDense = bfMatDenseRealToMatDense(matDenseReal);

  BfVecReal *rowView = bfVecRealNew();
  HANDLE_ERROR();

  BfSize n = bfMatDenseRealGetNumCols(matDenseReal);
  BfSize stride = bfMatDenseGetColStride(matDense);
  BfReal *data = matDenseReal->data + i*bfMatDenseGetRowStride(matDense);

  bfVecRealInitView(rowView, n, stride, data);

  BF_ERROR_END() {
    BF_DIE();
  }

  return rowView;
}

BfVec *bfMatDenseRealGetColView(BfMat *mat, BfSize j) {
  BF_ERROR_BEGIN();

  BfMatDense *matDense = NULL;
  BfMatDenseReal *matDenseReal = NULL;
  BfVecReal *colView = NULL;

  matDense = bfMatToMatDense(mat);
  HANDLE_ERROR();

  matDenseReal = bfMatToMatDenseReal(mat);
  HANDLE_ERROR();

  colView = bfVecRealNew();
  HANDLE_ERROR();

  BfSize m = bfMatGetNumRows(mat);
  BfSize stride = bfMatDenseGetRowStride(matDense);
  BfReal *data = matDenseReal->data + j*bfMatDenseGetColStride(matDense);

  bfVecRealInitView(colView, m, stride, data);

  BF_ERROR_END()
    bfVecRealDeinitAndDealloc(&colView);

  return bfVecRealToVec(colView);
}

void bfMatDenseRealDelete(BfMatDenseReal **matDenseReal) {
  bfMatDenseRealDeinitAndDealloc(matDenseReal);
}

BfMat *bfMatDenseRealZerosLike(BfMatDenseReal const *matDenseReal, BfSize numRows, BfSize numCols) {
  BF_ERROR_BEGIN();

  BfMatDenseReal *zeros = NULL;

  if (numRows == BF_SIZE_BAD_VALUE)
    numRows = bfMatDenseRealGetNumRows(matDenseReal);

  if (numCols == BF_SIZE_BAD_VALUE)
    numCols = bfMatDenseRealGetNumCols(matDenseReal);

  zeros = bfMatDenseRealNewZeros(numRows, numCols);
  HANDLE_ERROR();

  BF_ERROR_END()
    bfMatDenseRealDeinitAndDealloc(&zeros);

  return bfMatDenseRealToMat(zeros);
}

BfType bfMatDenseRealGetType(BfMat const *mat) {
  (void)mat;
  return BF_TYPE_MAT_DENSE_REAL;
}

BfSize bfMatDenseRealNumBytes(BfMatDenseReal const *matDenseReal) {
  BfMat const *mat = bfMatDenseRealConstToMatConst(matDenseReal);
  BfSize numRows = bfMatGetNumRows(mat);
  BfSize numCols = bfMatGetNumCols(mat);
  return sizeof(BfReal)*numRows*numCols;
}

void bfMatDenseRealSave(BfMatDenseReal const *matDenseReal, char const *path) {
  BF_ERROR_BEGIN();

  FILE *fp = fopen(path, "w");
  if (fp == NULL)
    RAISE_ERROR(BF_ERROR_FILE_ERROR);

  BfSize m = bfMatDenseRealGetNumRows(matDenseReal);
  BfSize n = bfMatDenseRealGetNumCols(matDenseReal);
  BfSize numElts = m*n;

  fwrite(matDenseReal->data, sizeof(BfReal), numElts, fp);
  /* TODO: error-handling */

  BF_ERROR_END() {}

  fclose(fp);
}

void bfMatDenseRealDump(BfMatDenseReal const *matDenseReal, FILE *fp) {
//   BF_ERROR_BEGIN();

  BfMat const *mat = bfMatDenseRealConstToMatConst(matDenseReal);

  /* Write the number of rows: */
  BfSize numRows = bfMatGetNumRows(mat);
  fwrite(&numRows, sizeof(BfSize), 1, fp);

  /* Write the number of columns: */
  BfSize numCols = bfMatGetNumCols(mat);
  fwrite(&numCols, sizeof(BfSize), 1, fp);

  /* Write the contents of the matrix: */
  for (BfSize i = 0; i < numRows; ++i) {
    BfReal const *readPtr = matDenseReal->data + i*matDenseReal->super.rowStride;
    for (BfSize j = 0; j < numCols; ++j) {
      fwrite(readPtr, sizeof(BfReal), 1, fp);
      readPtr += matDenseReal->super.colStride;
    }
  }

//   BF_ERROR_END() {
//     BF_DIE();
//   }
}

void bfMatDenseRealPrint(BfMat const *mat, FILE *fp) {
  BF_ERROR_BEGIN();

  BfMatDense const *matDense = NULL;
  BfMatDenseReal const *matDenseReal = NULL;

  if (fp == NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  matDense = bfMatConstToMatDenseConst(mat);
  HANDLE_ERROR();

  matDenseReal = bfMatConstToMatDenseRealConst(mat);
  HANDLE_ERROR();

  BfSize numRows = bfMatGetNumRows(mat);
  BfSize numCols = bfMatGetNumCols(mat);
  BfSize rowStride = bfMatDenseGetRowStride(matDense);
  BfSize colStride = bfMatDenseGetColStride(matDense);

  fprintf(fp, "[");
  for (BfSize i = 0; i < numRows; ++i) {
    BfReal const *rowPtr = matDenseReal->data + i*rowStride;
    fprintf(fp, "[");
    for (BfSize j = 0; j < numCols - 1; ++j) {
      fprintf(fp, "%g, ", *rowPtr);
      rowPtr += colStride;
    }
    fprintf(fp, "%g]", *rowPtr);
  }
  fprintf(fp, "]\n");

  BF_ERROR_END() {}
}

BfSize bfMatDenseRealGetNumRows(BfMatDenseReal const *matDenseReal) {
  BfMat const *mat = bfMatDenseRealConstToMatConst(matDenseReal);
  return bfMatIsTransposed(mat) ? mat->numCols : mat->numRows;
}

BfSize bfMatDenseRealGetNumCols(BfMatDenseReal const *matDenseReal) {
  BfMat const *mat = bfMatDenseRealConstToMatConst(matDenseReal);
  return bfMatIsTransposed(mat) ? mat->numRows : mat->numCols;
}

void bfMatDenseRealSetRow(BfMat *mat, BfSize i, BfVec const *row) {
  BF_ERROR_BEGIN();

  BfVecReal const *rowReal = NULL;
  BfMatDense *matDense = NULL;
  BfMatDenseReal *matDenseReal = NULL;

  BfSize m = bfMatGetNumRows(mat);
  BfSize n = bfMatGetNumCols(mat);

  if (i >= m)
    RAISE_ERROR(BF_ERROR_OUT_OF_RANGE);

  rowReal = bfVecConstToVecRealConst(row);
  HANDLE_ERROR();

  if (row->size > n)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  matDense = bfMatToMatDense(mat);
  HANDLE_ERROR();

  matDenseReal = bfMatToMatDenseReal(mat);
  HANDLE_ERROR();

  BfSize rowStride = bfMatDenseGetRowStride(matDense);
  BfSize colStride = bfMatDenseGetColStride(matDense);

  if (matDenseReal->data == NULL)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  BfReal *dstPtr = matDenseReal->data + i*rowStride;
  if (dstPtr == NULL)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  BfReal const *srcPtr = rowReal->data;
  if (srcPtr == NULL)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  for (BfSize j = 0; j < n; ++j) {
    *dstPtr = *srcPtr;
    dstPtr += colStride;
    srcPtr += rowReal->stride;
  }

  BF_ERROR_END() {}
}

void bfMatDenseRealSetCol(BfMatDenseReal *matDenseReal, BfSize j, BfVec const *col) {
  BF_ERROR_BEGIN();

  BfVecReal const *colReal = NULL;
  BfMatDense *matDense = NULL;

  BfSize m = bfMatDenseRealGetNumRows(matDenseReal);
  BfSize n = bfMatDenseRealGetNumCols(matDenseReal);

  if (j >= n)
    RAISE_ERROR(BF_ERROR_OUT_OF_RANGE);

  colReal = bfVecConstToVecRealConst(col);
  HANDLE_ERROR();

  if (col->size > m)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  matDense = bfMatDenseRealToMatDense(matDenseReal);
  HANDLE_ERROR();

  BfSize rowStride = bfMatDenseGetRowStride(matDense);
  BfSize colStride = bfMatDenseGetColStride(matDense);

  if (matDenseReal->data == NULL)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  BfReal *dstPtr = matDenseReal->data + j*colStride;
  if (dstPtr == NULL)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  BfReal const *srcPtr = colReal->data;
  if (srcPtr == NULL)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  for (BfSize i = 0; i < m; ++i) {
    *dstPtr = *srcPtr;
    dstPtr += rowStride;
    srcPtr += colReal->stride;
  }

  BF_ERROR_END() {}
}

BfMat *bfMatDenseRealGetRowRange(BfMat *mat, BfSize i0, BfSize i1) {
  BfSize numRows = mat->numRows;

  BF_ASSERT(i0 < i1);
  BF_ASSERT(i1 <= numRows);
  BF_ASSERT(!bfMatIsTransposed(mat)); // TODO: implement

  BF_ERROR_BEGIN();

  BfMat *matView = bfMatGetView(mat);

  BfMatDenseReal *submat = bfMatToMatDenseReal(matView);
  HANDLE_ERROR();

  if (i1 - i0 != numRows) {
    bfMatDenseRealToMat(submat)->numRows = i1 - i0;
    submat->data += bfMatDenseRealToMatDense(submat)->rowStride*i0;
  }

  BF_ERROR_END()
    bfMatDelete(&matView);

  return matView;
}

BfMat *bfMatDenseRealGetRowRangeCopy(BfMatDenseReal const *matDenseReal, BfSize i0, BfSize i1) {
  BF_ERROR_BEGIN();

  BfMat const *mat = bfMatDenseRealConstToMatConst(matDenseReal);
  HANDLE_ERROR();

  BfMatDense const *matDense = bfMatDenseRealConstToMatDenseConst(matDenseReal);
  HANDLE_ERROR();

  BfSize m = bfMatGetNumRows(mat);
  BfSize n = bfMatGetNumCols(mat);

  if (i0 >= i1)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  if (i1 > m)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  BfMatDenseReal *matDenseRealCopy = bfMatDenseRealNew();
  HANDLE_ERROR();

  bfMatDenseRealInit(matDenseRealCopy, i1 - i0, n);
  HANDLE_ERROR();

  BfMatDense *matDenseCopy = bfMatDenseRealToMatDense(matDenseRealCopy);
  HANDLE_ERROR();

  for (BfSize i = i0; i < i1; ++i) {
    BfReal const *readPtr = matDenseReal->data + i*matDense->rowStride;
    BfReal *writePtr = matDenseRealCopy->data + (i - i0)*matDenseCopy->rowStride;
    for (BfSize j = 0; j < n; ++j) {
      *writePtr = *readPtr;
      readPtr += matDense->colStride;
      writePtr += matDenseCopy->colStride;
    }
  }

  BF_ERROR_END() {
    bfMatDenseRealDeinitAndDealloc(&matDenseRealCopy);
  }

  return bfMatDenseRealToMat(matDenseRealCopy);
}

BfMat *bfMatDenseRealGetColRangeCopy(BfMatDenseReal const *matDenseReal, BfSize j0, BfSize j1) {
  BF_ERROR_BEGIN();

  if (j0 > j1)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfMat const *mat = bfMatDenseRealConstToMatConst(matDenseReal);

  BfMatDense const *matDense = bfMatDenseRealConstToMatDenseConst(matDenseReal);
  HANDLE_ERROR();

  BfSize m = bfMatGetNumRows(mat);
  BfSize n = bfMatGetNumCols(mat);

  if (j1 > n)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  BfMatDenseReal *matDenseRealCopy = bfMatDenseRealNew();
  HANDLE_ERROR();

  bfMatDenseRealInit(matDenseRealCopy, m, j1 - j0);
  HANDLE_ERROR();

  BfMatDense *matDenseCopy = bfMatDenseRealToMatDense(matDenseRealCopy);
  HANDLE_ERROR();

  for (BfSize j = j0; j < j1; ++j) {
    BfReal const *readPtr = matDenseReal->data + j*matDense->colStride;
    BfReal *writePtr = matDenseRealCopy->data + (j - j0)*matDenseCopy->colStride;
    for (BfSize i = 0; i < m; ++i) {
      *writePtr = *readPtr;
      readPtr += matDense->rowStride;
      writePtr += matDenseCopy->rowStride;
    }
  }

  BF_ERROR_END() {
    bfMatDenseRealDeinitAndDealloc(&matDenseRealCopy);
  }

  return bfMatDenseRealToMat(matDenseRealCopy);
}

void bfMatDenseRealPermuteRows(BfMatDenseReal *matDenseReal, BfPerm const *perm) {
  BF_ERROR_BEGIN();

  BfMat *mat = bfMatDenseRealToMat(matDenseReal);

  BfSize m = bfMatGetNumRows(mat);
  BfSize n = bfMatGetNumCols(mat);

  if (bfPermGetSize(perm) != m)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfMatDense *matDense = NULL;
  BfMatDense *matDensePerm = NULL;
  BfMatDenseReal *matDenseRealPerm = NULL;

  BfMat *matPerm = bfMatCopy(mat);
  HANDLE_ERROR();

  matDense = bfMatToMatDense(mat);
  HANDLE_ERROR();

  matDensePerm = bfMatToMatDense(matPerm);
  HANDLE_ERROR();

  matDenseReal = bfMatToMatDenseReal(mat);
  HANDLE_ERROR();

  matDenseRealPerm = bfMatToMatDenseReal(matPerm);
  HANDLE_ERROR();

  for (BfSize i = 0; i < m; ++i) {
    BfReal const *inRowPtr = matDenseRealPerm->data + i*matDensePerm->rowStride;
    BfReal *outRowPtr = matDenseReal->data + perm->index[i]*matDense->rowStride;
    for (BfSize j = 0; j < n; ++j) {
      *outRowPtr = *inRowPtr;
      inRowPtr += matDensePerm->colStride;
      outRowPtr += matDense->colStride;
    }
  }

  BF_ERROR_END() {
    BF_DIE();
  }

  bfMatDelete(&matPerm);
}

void bfMatDenseRealPermuteCols(BfMatDenseReal *matDenseReal, BfPerm const *perm) {
  BF_ERROR_BEGIN();

  BfMat *mat = bfMatDenseRealToMat(matDenseReal);

  BfSize m = bfMatGetNumRows(mat);
  BfSize n = bfMatGetNumCols(mat);

  if (bfPermGetSize(perm) != n)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfMatDense *matDense = NULL;
  BfMatDense *matDensePerm = NULL;
  BfMatDenseReal *matDenseRealPerm = NULL;

  BfMat *matPerm = bfMatCopy(mat);
  HANDLE_ERROR();

  matDense = bfMatToMatDense(mat);
  HANDLE_ERROR();

  matDensePerm = bfMatToMatDense(matPerm);
  HANDLE_ERROR();

  matDenseReal = bfMatToMatDenseReal(mat);
  HANDLE_ERROR();

  matDenseRealPerm = bfMatToMatDenseReal(matPerm);
  HANDLE_ERROR();

  for (BfSize j = 0; j < n; ++j) {
    BfReal const *inColPtr = matDenseRealPerm->data + j*matDensePerm->colStride;
    BfReal *outColPtr = matDenseReal->data + perm->index[j]*matDense->colStride;
    for (BfSize i = 0; i < m; ++i) {
      *outColPtr = *inColPtr;
      inColPtr += matDensePerm->rowStride;
      outColPtr += matDense->rowStride;
    }
  }

  BF_ERROR_END() {
    BF_DIE();
  }

  bfMatDelete(&matPerm);
}

/** Interface: MatDense */

static BfMatDenseVtable MAT_DENSE_VTABLE = {
  .Svd = (__typeof__(&bfMatDenseSvd))bfMatDenseRealSvd
};

/** Upcasting: */

BfMat *bfMatDenseRealToMat(BfMatDenseReal *matDenseReal) {
  return &matDenseReal->super.super;
}

BfMat const *bfMatDenseRealConstToMatConst(BfMatDenseReal const *matDenseReal) {
  return &matDenseReal->super.super;
}

BfMatDense *bfMatDenseRealToMatDense(BfMatDenseReal *matDenseReal) {
  return &matDenseReal->super;
}

BfMatDense const *bfMatDenseRealConstToMatDenseConst(BfMatDenseReal const *matDenseReal) {
  return &matDenseReal->super;
}

/** Downcasting: */

BfMatDenseReal *bfMatToMatDenseReal(BfMat *mat) {
  if (!bfMatInstanceOf(mat, BF_TYPE_MAT_DENSE_REAL)) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return NULL;
  } else {
    return (BfMatDenseReal *)mat;
  }
}

BfMatDenseReal const *bfMatConstToMatDenseRealConst(BfMat const *mat) {
  if (!bfMatInstanceOf(mat, BF_TYPE_MAT_DENSE_REAL)) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return NULL;
  } else {
    return (BfMatDenseReal const *)mat;
  }
}

/** Implementation: MatDenseReal */

BfMatDenseReal *bfMatDenseRealNew() {
  BF_ERROR_BEGIN();

  BfMatDenseReal *mat = bfMemAlloc(1, sizeof(BfMatDenseReal));
  if (mat == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BF_ERROR_END() {}

  return mat;
}

BfMatDenseReal *bfMatDenseRealNewViewFromPtr(BfSize numRows, BfSize numCols, BfReal *data, BfSize rowStride, BfSize colStride) {
  BF_ERROR_BEGIN();

  BfMatDenseReal *matDenseReal = bfMatDenseRealNew();
  HANDLE_ERROR();

  bfMatDenseRealInitViewFromPtr(matDenseReal, numRows, numCols, data, rowStride, colStride);

  BF_ERROR_END() {
    BF_DIE();
  }

  return matDenseReal;
}

BfMatDenseReal const *bfMatDenseRealNewViewConstFromPtr(BfSize numRows, BfSize numCols, BfReal const *data, BfSize rowStride, BfSize colStride) {
  BF_ERROR_BEGIN();

  BfMatDenseReal *matDenseReal = bfMatDenseRealNew();
  HANDLE_ERROR();

  bfMatDenseRealInitViewConstFromPtr(matDenseReal, numRows, numCols, data, rowStride, colStride);

  BF_ERROR_END() {
    BF_DIE();
  }

  return matDenseReal;
}

BfMatDenseReal *bfMatDenseRealNewWithValue(BfSize numRows, BfSize numCols, BfReal value) {
  BF_ERROR_BEGIN();

  BfMatDenseReal *matDenseReal = bfMemAlloc(1, sizeof(BfMatDenseReal));
  if (matDenseReal == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  bfMatDenseRealInit(matDenseReal, numRows, numCols);
  HANDLE_ERROR();

  for (BfSize i = 0; i < numRows; ++i) {
    BfReal *writePtr = matDenseReal->data + i*matDenseReal->super.rowStride;
    for (BfSize j = 0; j < numCols; ++j) {
      *writePtr = value;
      writePtr += matDenseReal->super.colStride;
    }
  }

  BF_ERROR_END() {
    BF_DIE();
  }

  return matDenseReal;
}

BfMatDenseReal *bfMatDenseRealNewFromMatrix(BfMat const *mat) {
  BF_ERROR_BEGIN();

  BfMatDenseReal *matDenseReal = bfMatDenseRealNew();
  HANDLE_ERROR();

  BfSize numRows = bfMatGetNumRows(mat);
  BfSize numCols = bfMatGetNumCols(mat);

  BfType matType = bfMatGetType(mat);

  /* If `mat` is some kind of block matrix, then a reasonable thing to
   * do is iterate over each of its blocks and set the corresponding
   * blocks of `matDenseReal` one at a time.
   *
   * NOTE: we can do a little better here by specializing on the type
   * of block matrix, but this is a reasonble and simple way to get
   * this conversion working. */
  if (bfTypeDerivedFrom(matType, BF_TYPE_MAT_BLOCK)) {
    bfMatDenseRealInit(matDenseReal, numRows, numCols);
    HANDLE_ERROR();

    BfMatBlock const *matBlock = bfMatConstToMatBlockConst(mat);

    BfSize numRowBlocks = bfMatBlockGetNumRowBlocks(matBlock);
    BfSize numColBlocks = bfMatBlockGetNumColBlocks(matBlock);

    for (BfSize i = 0; i < numRowBlocks; ++i) {
      BfSize i0 = bfMatBlockGetRowOffset(matBlock, i);
      BfSize i1 = i0 + bfMatBlockGetNumBlockRows(matBlock, i);

      for (BfSize j = 0; j < numColBlocks; ++j) {
        BfSize j0 = bfMatBlockGetColOffset(matBlock, j);
        BfSize j1 = j0 + bfMatBlockGetNumBlockCols(matBlock, j);

        BfMat const *block = bfMatBlockGetBlockConst(matBlock, i, j);

        BfMatDenseReal *block_ = bfMatDenseRealNewFromMatrix(block);
        HANDLE_ERROR();

        bfMatDenseRealSetBlock(matDenseReal, i0, i1, j0, j1, block_);

        bfMatDenseRealDeinitAndDealloc(&block_);

        bfMatDelete((BfMat **)&block);
      }
    }
  }

  else if (matType == BF_TYPE_MAT_DENSE_REAL) {
    bfMatDenseRealInit(matDenseReal, numRows, numCols);
    HANDLE_ERROR();

    BfMatDenseReal const *matDenseReal_ = bfMatConstToMatDenseRealConst(mat);
    bfMatDenseRealSetBlock(matDenseReal, 0, numRows, 0, numCols, matDenseReal_);
  }

  else if (matType == BF_TYPE_MAT_ZERO) {
    bfMatDenseRealInitWithValue(matDenseReal, numRows, numCols, 0.0);
    HANDLE_ERROR();
  }

  else if (matType == BF_TYPE_MAT_IDENTITY) {
    if (numRows != numCols)
      RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

    bfMatDenseRealInitWithValue(matDenseReal, numRows, numCols, 0.0);
    HANDLE_ERROR();

    for (BfSize i = 0; i < numRows; ++i) {
      *(matDenseReal->data + i*matDenseReal->super.rowStride
        + i*matDenseReal->super.colStride) = 1.0;
    }
  }

  /* TODO: Haven't implemented the conversion to `MatDenseReal` for
   * this type of matrix yet... */
  else {
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);
  }

  BF_ERROR_END() {
    bfMatDenseRealDeinitAndDealloc(&matDenseReal);

    BF_DIE();
  }

  return matDenseReal;
}

BfMatDenseReal *bfMatDenseRealNewZeros(BfSize numRows, BfSize numCols) {
  BF_ERROR_BEGIN();

  BfMatDenseReal *matDenseReal = bfMatDenseRealNew();
  HANDLE_ERROR();

  bfMatDenseRealInitZeros(matDenseReal, numRows, numCols);

  BF_ERROR_END() {
    BF_DIE();
  }

  return matDenseReal;
}

BfMatDenseReal *bfMatDenseRealFromFile(char const *path, BfSize numRows, BfSize numCols) {
  BF_ERROR_BEGIN();

  FILE *fp = NULL;
  BfMatDenseReal *matDenseReal = NULL;

  if (numRows == BF_SIZE_BAD_VALUE && numCols == BF_SIZE_BAD_VALUE)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  matDenseReal = bfMatDenseRealNew();
  HANDLE_ERROR();

  fp = fopen(path, "r");
  if (fp == NULL)
    RAISE_ERROR(BF_ERROR_FILE_ERROR);

  /* If we didn't pass the number of rows or columns, figure out the
   * matrix shape from the size of the binary file. */
  if (numRows == BF_SIZE_BAD_VALUE || numCols == BF_SIZE_BAD_VALUE) {
    BF_ASSERT(numRows != BF_SIZE_BAD_VALUE || numCols != BF_SIZE_BAD_VALUE);
    fseek(fp, 0, SEEK_END);
    BfSize numBytes = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    if (numBytes % sizeof(BfReal) != 0)
      RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);
    BfSize numElts = numBytes/sizeof(BfReal);
    if (numRows == BF_SIZE_BAD_VALUE) numRows = numElts/numCols;
    if (numCols == BF_SIZE_BAD_VALUE) numCols = numElts/numRows;
  }

  bfMatDenseRealInit(matDenseReal, numRows, numCols);
  HANDLE_ERROR();

  BfSize size = numRows*numCols;
  if (fread(matDenseReal->data, sizeof(BfReal), size, fp) != size)
    RAISE_ERROR(BF_ERROR_FILE_ERROR);

  BF_ERROR_END() {
    BF_DIE();
  }

  fclose(fp);

  return matDenseReal;
}

BfMatDenseReal *bfMatDenseRealNewFromRealArray(BfRealArray *realArray, BfSize numRows, BfSize numCols, BfPolicy policy) {
  BF_ERROR_BEGIN();

  BfMatDenseReal *matDenseReal = bfMatDenseRealNew();
  HANDLE_ERROR();

  bfMatDenseRealInitFromRealArray(matDenseReal, realArray, numRows, numCols, policy);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  return matDenseReal;
}

BfMatDenseReal *bfMatDenseRealNewFromCsv(char const *path) {
   BF_ERROR_BEGIN();

  BfMatDenseReal *matDenseReal = NULL;

  FILE *fp = NULL;

  char *lineptr = NULL;
  ssize_t nread;
  size_t n = 1024;
  char *saveptr;
  char *tok;

  BfRealArray *values = NULL;

  BfSize numRows = 0;
  BfSize numCols = BF_SIZE_BAD_VALUE;

  values = bfRealArrayNewWithDefaultCapacity();
  HANDLE_ERROR();

  fp = fopen(path, "r");
  if (fp == NULL)
    RAISE_ERROR(BF_ERROR_FILE_ERROR);

  /* Iterate over the lines in the CSV file, accumulating them into
  `values. If we discover that the lines have different numbers of
  values, we bail. */
  while (true) {
    lineptr = NULL;
    saveptr = NULL;

    nread = getline(&lineptr, &n, fp);
    if (nread < 0)
      break;

    BfSize lineNumCols = 0;
    while (true) {
      tok = strtok_r(lineNumCols ? NULL : lineptr, " ", &saveptr);
      if (tok == NULL)
        break;

      errno = 0;
      double value = strtod(tok, NULL);
      if (errno != 0)
        RAISE_ERROR(BF_ERROR_FILE_ERROR);

      bfRealArrayAppend(values, value);
      HANDLE_ERROR();

      ++lineNumCols;
    }

    if (!BF_SIZE_OK(numCols)) numCols = lineNumCols;
    if (lineNumCols != numCols)
      RAISE_ERROR(BF_ERROR_FILE_ERROR);

    ++numRows;
  }

  bfRealArrayShrinkCapacityToSize(values);
  HANDLE_ERROR();

  matDenseReal = bfMatDenseRealNewFromRealArray(values, numRows, numCols, BF_POLICY_STEAL);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  fclose(fp);

  bfRealArrayDeinitAndDealloc(&values);

  return matDenseReal;
}

void bfMatDenseRealInit(BfMatDenseReal *matDenseReal, BfSize numRows, BfSize numCols) {
  BF_ERROR_BEGIN();

  bfMatDenseInit(bfMatDenseRealToMatDense(matDenseReal),
                 &MAT_VTABLE, &MAT_DENSE_VTABLE, numRows, numCols, numCols, 1);

  matDenseReal->data = bfMemAlloc(numRows*numCols, sizeof(BfReal));
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }
}

void bfMatDenseRealInitViewFromPtr(BfMatDenseReal *matDenseReal, BfSize numRows, BfSize numCols, BfReal *data, BfSize rowStride, BfSize colStride) {
  BF_ERROR_BEGIN();

  bfMatDenseInit(&matDenseReal->super, &MAT_VTABLE, &MAT_DENSE_VTABLE, numRows, numCols, rowStride, colStride);
  HANDLE_ERROR();

  matDenseReal->super.super.props |= BF_MAT_PROPS_VIEW;

  matDenseReal->data = data;

  BF_ERROR_END() {
    BF_DIE();
  }
}

void bfMatDenseRealInitViewConstFromPtr(BfMatDenseReal const *matDenseReal, BfSize numRows, BfSize numCols, BfReal const *data, BfSize rowStride, BfSize colStride) {
  BF_ERROR_BEGIN();

  /* Get writable pointer and use that: */
  BfMatDenseReal *matDenseReal_ = (BfMatDenseReal *)matDenseReal;

  bfMatDenseInit(&matDenseReal_->super, &MAT_VTABLE, &MAT_DENSE_VTABLE, numRows, numCols, rowStride, colStride);
  HANDLE_ERROR();

  matDenseReal_->super.super.props |= BF_MAT_PROPS_VIEW;

  matDenseReal_->data = (BfReal *)data;

  BF_ERROR_END() {
    BF_DIE();
  }
}

void bfMatDenseRealInitCopy(BfMatDenseReal *matDenseReal,
                            BfMatDenseReal const *otherMatDenseReal) {
  BF_ERROR_BEGIN();

  BfMat const *otherMat = bfMatDenseRealConstToMatConst(otherMatDenseReal);

  BfSize m = bfMatGetNumRows(otherMat);
  BfSize n = bfMatGetNumCols(otherMat);

  bfMatDenseRealInit(matDenseReal, m, n);
  HANDLE_ERROR();

  BfMatDense *matDense = bfMatDenseRealToMatDense(matDenseReal);
  BfSize rowStride = bfMatDenseGetRowStride(matDense);
  BfSize colStride = bfMatDenseGetColStride(matDense);

  BfMatDense const *otherMatDense = bfMatDenseRealConstToMatDenseConst(otherMatDenseReal);
  BfSize otherRowStride = bfMatDenseGetRowStride(otherMatDense);
  BfSize otherColStride = bfMatDenseGetColStride(otherMatDense);

  for (BfSize i = 0; i < m; ++i) {
    BfReal *outPtr = matDenseReal->data + i*rowStride;
    BfReal const *inPtr = otherMatDenseReal->data + i*otherRowStride;
    for (BfSize j = 0; j < n; ++j) {
      *outPtr = *inPtr;
      outPtr += colStride;
      inPtr += otherColStride;
    }
  }

  BF_ERROR_END() {}
}

void bfMatDenseRealInitFromPtr(BfMatDenseReal *matDenseReal, BfSize numRows, BfSize numCols, BfReal *data, BfPolicy policy) {
  BF_ERROR_BEGIN();

  bfMatDenseInit(bfMatDenseRealToMatDense(matDenseReal),
                 &MAT_VTABLE, &MAT_DENSE_VTABLE, numRows, numCols, numCols, 1);

  if (policy == BF_POLICY_COPY) {
    matDenseReal->data = bfMemAlloc(numRows*numCols, sizeof(BfReal));
    HANDLE_ERROR();

    bfMemCopy(data, numRows*numCols, sizeof(BfReal), matDenseReal->data);
  }

  else if (policy == BF_POLICY_VIEW || policy == BF_POLICY_STEAL) {
    matDenseReal->data = data;
  }

  else BF_DIE();

  BF_ERROR_END() {
    BF_DIE();
  }
}

void bfMatDenseRealInitFromRealArray(BfMatDenseReal *matDenseReal, BfRealArray *realArray, BfSize numRows, BfSize numCols, BfPolicy policy) {
  BF_ERROR_BEGIN();

  if (bfRealArrayGetSize(realArray) != numRows*numCols)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (policy == BF_POLICY_STEAL) {
    if (realArray->isView)
      RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);
    realArray->isView = true;
  }

  bfMatDenseRealInitFromPtr(matDenseReal, numRows, numCols, realArray->data, BF_POLICY_STEAL);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }
}

void bfMatDenseRealInitWithValue(BfMatDenseReal *matDenseReal, BfSize numRows,
                                 BfSize numCols, BfReal fillValue) {
  BF_ERROR_BEGIN();

  BfMatDense *matDense = NULL;

  bfMatDenseRealInit(matDenseReal, numRows, numCols);
  HANDLE_ERROR();

  matDense = bfMatDenseRealToMatDense(matDenseReal);
  HANDLE_ERROR();

  BfSize colStride = bfMatDenseGetColStride(matDense);

  BfReal *ptr = matDenseReal->data;
  for (BfSize i = 0; i < numRows; ++i) {
    for (BfSize j = 0; j < numCols; ++j) {
      *ptr = fillValue;
      ptr += colStride;
    }
  }

  BF_ERROR_END()
    bfMatDenseRealDeinit(matDenseReal);
}

void bfMatDenseRealInitZeros(BfMatDenseReal *matDenseReal, BfSize numRows, BfSize numCols) {
  BF_ERROR_BEGIN();

  bfMatDenseRealInit(matDenseReal, numRows, numCols);
  HANDLE_ERROR();

  memset(matDenseReal->data, 0x0, numRows*numCols*sizeof(BfReal));

  BF_ERROR_END() {
    BF_DIE();
  }
}

void bfMatDenseRealDeinit(BfMatDenseReal *matDenseReal) {
  BfMat *mat = bfMatDenseRealToMat(matDenseReal);

  if (!bfMatIsView(mat))
    free(matDenseReal->data);

  matDenseReal->data = NULL;

  bfMatDenseDeinit(&matDenseReal->super);
}

void bfMatDenseRealDealloc(BfMatDenseReal **mat) {
  free(*mat);
  *mat = NULL;
}

void bfMatDenseRealDeinitAndDealloc(BfMatDenseReal **mat) {
  bfMatDenseRealDeinit(*mat);
  bfMatDenseRealDealloc(mat);
}

void bfMatDenseRealSvd(BfMatDenseReal const *mat, BfMatDenseReal **UPtr,
                       BfMatDiagReal **SPtr, BfMatDenseReal **VTPtr) {
  BF_ERROR_BEGIN();

  BfMat const *super = bfMatDenseRealConstToMatConst(mat);

  BfSize m = super->numRows;
  BfSize n = super->numCols;
  BfSize k = m < n ? m : n;

  BfMatDenseReal *U = bfMatDenseRealNew();
  HANDLE_ERROR();

  bfMatDenseRealInit(U, m, k);
  HANDLE_ERROR();

  BfMatDiagReal *S = bfMatDiagRealNew();
  HANDLE_ERROR();

  bfMatDiagRealInit(S, k, k);
  HANDLE_ERROR();

  BfMatDenseReal *VT = bfMatDenseRealNew();
  HANDLE_ERROR();

  bfMatDenseRealInit(VT, k, n);
  HANDLE_ERROR();

  BfReal *superb = NULL;

  /* dgesvd will overwrite A, so allocate space for a copy */
  BfReal *dataCopy = bfMemAlloc(m*n, sizeof(BfReal));
  if (dataCopy == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  /* copy contents of A */
  bfMemCopy(mat->data, m*n, sizeof(BfReal), dataCopy);

  /* output array which contains information about superdiagonal
   * elements which didn't converge
   *
   * more info here: tinyurl.com/2p8f5ev3 */
  superb = bfMemAlloc(((m < n) ? m : n) - 1, sizeof(BfReal));
  if (superb == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  /* compute the SVD */
  lapack_int info = LAPACKE_dgesvd(
    LAPACK_ROW_MAJOR, 'S', 'S', m, n, dataCopy, n, S->data,
    U->data, m < n ? m : n, VT->data, n, superb);

  /* check for invalid arguments */
  if (info < 0)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  /* check for errors */
  if (info > 0)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  /* Make these both orthogonal: */
  bfMatDenseRealToMat(U)->props |= BF_MAT_PROPS_ORTHO;
  bfMatDenseRealToMat(VT)->props |= BF_MAT_PROPS_ORTHO;

  /* Update out pointers: */
  *UPtr = U;
  *SPtr = S;
  *VTPtr = VT;

  BF_ERROR_END() {
    BF_DIE();
  }

  free(dataCopy);
  free(superb);
}

void bfMatDenseRealSetBlock(BfMatDenseReal *matDenseReal,
                            BfSize i0, BfSize i1, BfSize j0, BfSize j1,
                            BfMatDenseReal const *block) {
  BF_ERROR_BEGIN();

  if (i0 > i1)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (j0 > j1)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  { BfMat *mat = bfMatDenseRealToMat(matDenseReal);

    if (bfMatGetNumRows(mat) < i1)
      RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

    if (bfMatGetNumCols(mat) < j1)
      RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS); }

  BfSize numBlockRows = i1 - i0;
  BfSize numBlockCols = j1 - j0;

  { BfMat const *mat = bfMatDenseRealConstToMatConst(block);

    if (bfMatGetNumRows(mat) != numBlockRows)
      RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

    if (bfMatGetNumCols(mat) != numBlockCols)
      RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS); }

  for (BfSize i = i0; i < i1; ++i) {
    BfReal *writePtr = matDenseReal->data
      + i*matDenseReal->super.rowStride + j0*matDenseReal->super.colStride;
    BfReal const *readPtr = block->data + (i - i0)*block->super.rowStride;
    for (BfSize j = 0; j < numBlockCols; ++j) {
      *writePtr = *readPtr;
      writePtr += matDenseReal->super.colStride;
      readPtr += block->super.colStride;
    }
  }

  BF_ERROR_END() {}
}

static void scaleRows_vecReal(BfMatDenseReal *matDenseReal, BfVecReal const *vecReal) {
  BF_ERROR_BEGIN();

  BfMat *mat = bfMatDenseRealToMat(matDenseReal);

  BfSize m = bfMatGetNumRows(mat);
  if (m != vecReal->super.size)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize n = bfMatGetNumCols(mat);

  for (BfSize i = 0; i < m; ++i) {
    BfReal *outPtr = matDenseReal->data + i*matDenseReal->super.rowStride;
    for (BfSize j = 0; j < n; ++j) {
      *outPtr *= vecReal->data[i];
      outPtr += matDenseReal->super.colStride;
    }
  }

  BF_ERROR_END() {}
}

void bfMatDenseRealScaleRows(BfMatDenseReal *matDenseReal, BfVec const *vec) {
  switch (bfVecGetType(vec)) {
  case BF_TYPE_VEC_REAL:
    scaleRows_vecReal(matDenseReal, bfVecConstToVecRealConst(vec));
    break;
  default:
    bfSetError(BF_ERROR_NOT_IMPLEMENTED);
    break;
  }
}

static BfMat *mul_matDiagReal(BfMatDenseReal const *matDenseReal,
                              BfMatDiagReal const *otherMatDiagReal) {
  BF_ERROR_BEGIN();

  BfMat const *mat = bfMatDenseRealConstToMatConst(matDenseReal);
  BfMat const *otherMat = bfMatDiagRealConstToMatConst(otherMatDiagReal);
  BfMatDenseReal *newMatDenseReal = NULL;

  BfSize m = bfMatGetNumRows(mat);
  BfSize n = bfMatGetNumCols(mat);
  BfSize p = bfMatGetNumCols(otherMat);

  if (n != bfMatGetNumRows(otherMat))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  newMatDenseReal = bfMatDenseRealNew();
  HANDLE_ERROR();

  bfMatDenseRealInit(newMatDenseReal, m, p);
  HANDLE_ERROR();

  for (BfSize i = 0; i < m; ++i) {
    BfReal const *inPtr = matDenseReal->data + i*matDenseReal->super.rowStride;
    BfReal *outPtr = newMatDenseReal->data + i*newMatDenseReal->super.rowStride;
    for (BfSize j = 0; j < p; ++j) {
      *outPtr = *inPtr * otherMatDiagReal->data[j];
      inPtr += matDenseReal->super.colStride;
      outPtr += newMatDenseReal->super.colStride;
    }
  }

  BF_ERROR_END() {
    BF_DIE();
  }

  return bfMatDenseRealToMat(newMatDenseReal);
}

BfMat *bfMatDenseRealMul(BfMatDenseReal const *matDenseReal, BfMat const *otherMat) {
  switch (bfMatGetType(otherMat)) {
  case BF_TYPE_MAT_DIAG_REAL:
    return mul_matDiagReal(matDenseReal, bfMatConstToMatDiagRealConst(otherMat));
  default:
    bfSetError(BF_ERROR_NOT_IMPLEMENTED);
    return NULL;
  }
}

static BfVec *mulVec_vecReal(BfMatDenseReal const *matDenseReal,
                             BfVecReal const *vecReal,
                             BfSize m, BfSize n) {
  BF_ERROR_BEGIN();

  enum CBLAS_TRANSPOSE trans = getCblasTranspose(matDenseReal);

  BfVecReal *result = bfVecRealNew();
  HANDLE_ERROR();

  bfVecRealInit(result, m);
  HANDLE_ERROR();

  BfSize lda = getLeadingDimension(matDenseReal);

  BfReal alpha = 1, beta = 0;

  if (lda > 0) {
    cblas_dgemv(CblasRowMajor, trans, m, n, alpha,
                matDenseReal->data, lda, vecReal->data,
                vecReal->stride, beta, result->data,
                result->stride);
  } else {
    BF_DIE();
  }

  BF_ERROR_END() {
    BF_DIE();
  }

  return bfVecRealToVec(result);
}

BfVec *bfMatDenseRealMulVec(BfMatDenseReal const *matDenseReal, BfVec const *vec) {
  BF_ERROR_BEGIN();

  BfMat const *mat = bfMatDenseRealConstToMatConst(matDenseReal);
  BfVec *result = NULL;

  BfSize m = bfMatGetNumRows(mat);
  BfSize n = bfMatGetNumCols(mat);

  if (n != vec->size)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (n == 0)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  matDenseReal = bfMatConstToMatDenseRealConst(mat);
  HANDLE_ERROR();

  switch (bfVecGetType(vec)) {
  case BF_TYPE_VEC_REAL:
    result = mulVec_vecReal(
      matDenseReal, bfVecConstToVecRealConst(vec), m, n);
    break;
  default:
    result = NULL;
    bfSetError(BF_ERROR_NOT_IMPLEMENTED);
    break;
  }

  BF_ERROR_END() {
    BF_DIE();
  }

  return result;
}

static BfVec *rmulVec_vecReal(BfMatDenseReal const *matDenseReal,
                              BfVecReal const *vecReal,
                              BfSize m, BfSize n) {
  BF_ERROR_BEGIN();

  enum CBLAS_TRANSPOSE trans =
    getCblasTranspose(matDenseReal) == CblasNoTrans ? CblasTrans : CblasNoTrans;

  BfVecReal *result = bfVecRealNew();
  HANDLE_ERROR();

  bfVecRealInit(result, n);
  HANDLE_ERROR();

  BfSize lda = getLeadingDimension(matDenseReal);

  BfReal alpha = 1, beta = 0;

  if (lda > 0) {
    cblas_dgemv(CblasRowMajor, trans, m, n, alpha,
                matDenseReal->data, lda, vecReal->data,
                vecReal->stride, beta, result->data,
                result->stride);
  } else {
    BF_DIE();
  }

  BF_ERROR_END() {
    BF_DIE();
  }

  return bfVecRealToVec(result);
}

// static void mulInplace_matDenseReal_matDenseReal(BfMatDenseReal const *op1, BfMatDenseReal const *op2, BfMatDenseReal *res) {
//   BF_ERROR_BEGIN();

//   if (bfMatIsTransposed(TO_MAT(res)))
//     RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);

//   enum CBLAS_TRANSPOSE transa = getCblasTranspose(op1);
//   enum CBLAS_TRANSPOSE transb = getCblasTranspose(op2);

//   BfSize m = bfMatDenseRealGetNumRows(res);
//   BfSize n = bfMatDenseRealGetNumCols(res);
//   BfSize k = bfMatDenseRealGetNumCols(op1);

//   BfReal alpha = 1;
//   BfReal const *op1Data = op1->data;
//   BfSize op1LeadingDim = getLeadingDimension(op1);

//   BfReal const *op2Data = op2->data;
//   BfSize op2LeadingDim = getLeadingDimension(op2);

//   BfReal beta = 0;
//   BfReal *resData = res->data;
//   BfSize resLeadingDim = getLeadingDimension(res);

//   cblas_dgemm(CblasRowMajor, transa, transb, m, n, k, alpha, op1Data, op1LeadingDim, op2Data, op2LeadingDim, beta, resData, resLeadingDim);

//   // TODO: finish implementing...
//   BF_DIE();

//   BF_ERROR_END() {
//     BF_DIE();
//   }
// }

void bfMatDenseRealMulInplace(BfMatDenseReal const *op1, BfMat const *op2, BfMat *res) {
  BF_ERROR_BEGIN();

  BfSize m = bfMatDenseRealGetNumRows(op1);
  BfSize n = bfMatDenseRealGetNumCols(op1);
  BfSize p = bfMatGetNumCols(op2);

  if (n != bfMatGetNumRows(op2))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (m != bfMatGetNumRows(res))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (p != bfMatGetNumCols(res))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfType op2Type = bfMatGetType(op2);
  BfType resType = bfMatGetType(res);
  if (op2Type == BF_TYPE_MAT_DENSE_REAL && resType == BF_TYPE_MAT_DENSE_REAL) {
    // mulInplace_matDenseReal_matDenseReal(op1, TO_MAT_DENSE_REAL(op2), TO_MAT_DENSE_REAL(res));
    BF_DIE();
    HANDLE_ERROR();
  } else {
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);
  }

  BF_ERROR_END() {
    BF_DIE();
  }
}

BfVec *bfMatDenseRealRmulVec(BfMatDenseReal const *matDenseReal, BfVec const *vec) {
  BF_ERROR_BEGIN();

  BfMat const *mat = bfMatDenseRealConstToMatConst(matDenseReal);
  BfVec *result = NULL;

  BfSize m = bfMatGetNumRows(mat);
  BfSize n = bfMatGetNumCols(mat);

  if (m != vec->size)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (m == 0)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  matDenseReal = bfMatConstToMatDenseRealConst(mat);
  HANDLE_ERROR();

  switch (bfVecGetType(vec)) {
  case BF_TYPE_VEC_REAL:
    result = rmulVec_vecReal(
      matDenseReal, bfVecConstToVecRealConst(vec), m, n);
    break;
  default:
    result = NULL;
    bfSetError(BF_ERROR_NOT_IMPLEMENTED);
    break;
  }

  BF_ERROR_END() {
    BF_DIE();
  }

  return result;
}

void bfMatDenseRealPrintBlocksDeep(BfMatDenseReal const *matDenseReal, FILE *fp, BfSize i0, BfSize j0, BfSize depth) {
  BfMat const *mat = bfMatDenseRealConstToMatConst(matDenseReal);

  BfSize i1 = i0 + bfMatGetNumRows(mat);
  BfSize j1 = j0 + bfMatGetNumCols(mat);

  fprintf(fp, "%u %lu %lu %lu %lu %lu\n", BF_TYPE_MAT_DENSE_REAL, i0, i1, j0, j1, depth);
}

BfReal bfMatDenseRealNormMax(BfMatDenseReal const *matDenseReal) {
  BfSize numRows = bfMatDenseRealGetNumRows(matDenseReal);
  BfSize numCols = bfMatDenseRealGetNumCols(matDenseReal);

  BfMatDense const *matDense = bfMatDenseRealConstToMatDenseConst(matDenseReal);
  BfSize rowStride = bfMatDenseGetRowStride(matDense);
  BfSize colStride = bfMatDenseGetColStride(matDense);

  BfReal norm = -BF_INFINITY;
  for (BfSize i = 0; i < numRows; ++i) {
    BfReal const *rowPtr = matDenseReal->data + i*rowStride;
    for (BfSize j = 0; j < numCols; ++j) {
      norm = fmax(fabs(*rowPtr), norm);
      rowPtr += colStride;
    }
  }

  return norm;
}

BfReal bfMatDenseRealDistMax(BfMatDenseReal const *matDenseReal, BfMatDenseReal const *otherMatDenseReal) {
  BF_ERROR_BEGIN();

  BfSize numRows = bfMatDenseRealGetNumRows(matDenseReal);
  BfSize otherNumRows = bfMatDenseRealGetNumRows(otherMatDenseReal);
  if (numRows != otherNumRows)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize numCols = bfMatDenseRealGetNumCols(matDenseReal);
  BfSize otherNumCols = bfMatDenseRealGetNumCols(otherMatDenseReal);
  if (numCols != otherNumCols)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfMatDense const *matDense = bfMatDenseRealConstToMatDenseConst(matDenseReal);
  BfSize rowStride = bfMatDenseGetRowStride(matDense);
  BfSize colStride = bfMatDenseGetColStride(matDense);

  BfMatDense const *otherMatDense = bfMatDenseRealConstToMatDenseConst(otherMatDenseReal);
  BfSize otherRowStride = bfMatDenseGetRowStride(otherMatDense);
  BfSize otherColStride = bfMatDenseGetColStride(otherMatDense);

  BfReal dist = -BF_INFINITY;
  for (BfSize i = 0; i < numRows; ++i) {
    BfReal const *rowPtr = matDenseReal->data + i*rowStride;
    BfReal const *otherRowPtr = otherMatDenseReal->data + i*otherRowStride;
    for (BfSize j = 0; j < numCols; ++j) {
      dist = fmax(fabs(*rowPtr - *otherRowPtr), dist);
      rowPtr += colStride;
      otherRowPtr += otherColStride;
    }
  }

  BF_ERROR_END() {
    BF_DIE();
  }

  return dist;
}
