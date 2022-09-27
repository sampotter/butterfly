#include <bf/mat_dense_real.h>

#include <assert.h>
#include <stdlib.h>

#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/vec_real.h>

/** Interface: Mat */

#define INTERFACE BF_INTERFACE_Mat
BF_DEFINE_VTABLE(Mat, MatDenseReal)
#undef INTERFACE

BF_STUB(BfMat *, MatDenseRealCopy, BfMat const *)
BF_STUB(BfMat *, MatDenseRealGetView, BfMat *)
BF_STUB(BfVec *, MatDenseRealGetRowView, BfMat *, BfSize)

BfVec *bfMatDenseRealGetColView(BfMat *mat, BfSize j) {
  BEGIN_ERROR_HANDLING();

  BfMatDenseReal *matDenseReal = NULL;
  BfVecReal *colView = NULL;

  matDenseReal = bfMatToMatDenseReal(mat);
  HANDLE_ERROR();

  colView = bfVecRealNew();
  HANDLE_ERROR();

  BfSize m = bfMatGetNumRows(mat);
  BfSize stride = matDenseReal->rowStride;
  BfReal *data = matDenseReal->data + j*matDenseReal->colStride;

  bfVecRealInitView(colView, m, stride, data);

  END_ERROR_HANDLING()
    bfVecRealDeinitAndDealloc(&colView);

  return bfVecRealToVec(colView);
}

BF_STUB(BfVec *, MatDenseRealGetColRangeView, BfMat *, BfSize, BfSize, BfSize)
BF_STUB(void, MatDenseRealDelete, BfMat **)
BF_STUB(BfMat *, MatDenseRealEmptyLike, BfMat const *, BfSize, BfSize)
BF_STUB(BfMat *, MatDenseRealZerosLike, BfMat const *, BfSize, BfSize)

BfType bfMatDenseRealGetType(BfMat const *mat) {
  (void)mat;
  return BF_TYPE_MAT_DENSE_REAL;
}

BF_STUB(bool, MatDenseRealInstanceOf, BfMat const *, BfType)
BF_STUB(BfSize, MatDenseRealNumBytes, BfMat const *)
BF_STUB(void, MatDenseRealSave, BfMat const *, char const *)

void bfMatDenseRealPrint(BfMat const *mat, FILE *fp) {
  BEGIN_ERROR_HANDLING();

  if (fp == NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfMatDenseReal const *matDenseReal = bfMatConstToMatDenseRealConst(mat);
  HANDLE_ERROR();

  BfSize numRows = bfMatGetNumRows(mat);
  BfSize numCols = bfMatGetNumCols(mat);

  fprintf(fp, "[");
  for (BfSize i = 0; i < numRows; ++i) {
    BfReal const *rowPtr = matDenseReal->data + i*matDenseReal->rowStride;
    fprintf(fp, "[");
    for (BfSize j = 0; j < numCols - 1; ++j) {
      fprintf(fp, "%g, ", *rowPtr);
      rowPtr += matDenseReal->colStride;
    }
    fprintf(fp, "%g]", *rowPtr);
  }
  fprintf(fp, "]\n");

  END_ERROR_HANDLING() {}
}

BfSize bfMatDenseRealGetNumRows(BfMat const *mat) {
  return bfMatIsTransposed(mat) ? mat->numCols : mat->numRows;
}

BfSize bfMatDenseRealGetNumCols(BfMat const *mat) {
  return bfMatIsTransposed(mat) ? mat->numRows : mat->numCols;
}

void bfMatDenseRealSetRow(BfMat *mat, BfSize i, BfVec const *row) {
  BEGIN_ERROR_HANDLING();

  BfSize m = bfMatGetNumRows(mat);
  BfSize n = bfMatGetNumCols(mat);

  if (i >= m)
    RAISE_ERROR(BF_ERROR_OUT_OF_RANGE);

  BfVecReal const *rowReal = bfVecConstToVecRealConst(row);
  HANDLE_ERROR();

  if (row->size > n)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfMatDenseReal *matDenseReal = bfMatToMatDenseReal(mat);
  HANDLE_ERROR();

  if (matDenseReal->data == NULL)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  BfReal *dstPtr = matDenseReal->data + i*matDenseReal->rowStride;
  if (dstPtr == NULL)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  BfReal const *srcPtr = rowReal->data;
  if (srcPtr == NULL)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  for (BfSize j = 0; j < n; ++j) {
    *dstPtr = *srcPtr;
    dstPtr += matDenseReal->colStride;
    srcPtr += rowReal->stride;
  }

  END_ERROR_HANDLING() {}
}

BF_STUB(void, MatDenseRealSetCol, BfMat *, BfSize, BfVec const *)
BF_STUB(void, MatDenseRealSetColRange, BfMat *, BfSize, BfSize, BfSize, BfVec const *)
BF_STUB(BfMat *, MatDenseRealGetRowRange, BfMat *, BfSize, BfSize)
BF_STUB(BfMat *, MatDenseRealGetColRange, BfMat *, BfSize, BfSize)
BF_STUB(BfMat *, MatDenseRealGetRowRangeCopy, BfMat const *, BfSize, BfSize)
BF_STUB(BfMat *, MatDenseRealGetColRangeCopy, BfMat const *, BfSize, BfSize)
BF_STUB(void, MatDenseRealSetRowRange, BfMat *, BfSize, BfSize, BfMat const *)
BF_STUB(void, MatDenseRealPermuteRows, BfMat *, BfPerm const *)
BF_STUB(void, MatDenseRealPermuteCols, BfMat *, BfPerm const *)
BF_STUB(BfVec *, MatDenseRealRowDists, BfMat const *, BfMat const *)
BF_STUB(BfVec *, MatDenseRealColDists, BfMat const *, BfMat const *)
BF_STUB(BfVec *, MatDenseRealColDots, BfMat const *, BfMat const *)
BF_STUB(BfVec *, MatDenseRealColNorms, BfMat const *)
BF_STUB(void, MatDenseRealScaleCols, BfMat *, BfVec const *)
BF_STUB(BfVec *, MatDenseRealSumCols, BfMat const *)
BF_STUB(void, MatDenseRealAddInplace, BfMat *, BfMat const *)
BF_STUB(void, MatDenseRealAddDiag, BfMat *, BfMat const *)
BF_STUB(BfMat *, MatDenseRealSub, BfMat const *, BfMat const *)
BF_STUB(void, MatDenseRealSubInplace, BfMat *, BfMat const *)
BF_STUB(BfMat *, MatDenseRealMul, BfMat const *, BfMat const *)
BF_STUB(BfVec *, MatDenseRealMulVec, BfMat const *, BfVec const *)
BF_STUB(void, MatDenseRealMulInplace, BfMat *, BfMat const *)
BF_STUB(BfMat *, MatDenseRealSolveLU, BfMat const *, BfMat const *)
BF_STUB(BfMat *, MatDenseRealLstSq, BfMat const *, BfMat const *)
BF_STUB(BfMat *, MatDenseRealGetGivensRotation, BfVec const *, BfSize, BfSize)
BF_STUB(bool, MatDenseRealIsUpperTri, BfMat const *)
BF_STUB(BfVec *, MatDenseRealBackwardSolveVec, BfMat const *, BfVec const *)

/** Upcasting: */

BfMat *bfMatDenseRealToMat(BfMatDenseReal *matDenseReal) {
  return &matDenseReal->super;
}

BfMat const *bfMatDenseRealConstToMatConst(BfMatDenseReal const *matDenseReal) {
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
  BEGIN_ERROR_HANDLING();

  BfMatDenseReal *mat = malloc(sizeof(BfMatDenseReal));
  if (mat == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING() {}

  return mat;
}

BfMatDenseReal *bfMatDenseRealFromFile(char const *path, BfSize numRows, BfSize numCols) {
  BEGIN_ERROR_HANDLING();

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
    assert(numRows != BF_SIZE_BAD_VALUE || numCols != BF_SIZE_BAD_VALUE);
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

  END_ERROR_HANDLING()
    bfMatDenseRealDeinitAndDealloc(&matDenseReal);

  return matDenseReal;
}

void bfMatDenseRealInit(BfMatDenseReal *mat, BfSize numRows, BfSize numCols) {
  BEGIN_ERROR_HANDLING();

  bfMatInit(&mat->super, &MatVtbl, numRows, numCols);
  HANDLE_ERROR();

  mat->rowStride = numCols;
  mat->colStride = 1;

  mat->data = malloc(numRows*numCols*sizeof(BfReal));
  if (mat->data == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING()
    bfMatDeinit(&mat->super);
}

void bfMatDenseRealInitWithValue(BfMatDenseReal *mat, BfSize numRows,
                                 BfSize numCols, BfReal fillValue) {
  BEGIN_ERROR_HANDLING();

  bfMatDenseRealInit(mat, numRows, numCols);
  HANDLE_ERROR();

  BfReal *ptr = mat->data;
  for (BfSize i = 0; i < numRows; ++i) {
    for (BfSize j = 0; j < numCols; ++j) {
      *ptr = fillValue;
      ptr += mat->colStride;
    }
  }

  END_ERROR_HANDLING()
    bfMatDenseRealDeinit(mat);
}

void bfMatDenseRealDeinit(BfMatDenseReal *mat) {
  if (!(mat->super.props & BF_MAT_PROPS_VIEW))
    free(mat->data);

  mat->data = NULL;
}

void bfMatDenseRealDealloc(BfMatDenseReal **mat) {
  free(*mat);
  *mat = NULL;
}

void bfMatDenseRealDeinitAndDealloc(BfMatDenseReal **mat) {
  bfMatDenseRealDeinit(*mat);
  bfMatDenseRealDealloc(mat);
}
