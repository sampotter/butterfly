#include <bf/mat_dense_real.h>

#include <assert.h>
#include <stdlib.h>

#include <bf/error.h>
#include <bf/error_macros.h>

/** Interface: Mat */

#define INTERFACE BF_INTERFACE_Mat
BF_DEFINE_VTABLE(Mat, MatDenseReal)
#undef INTERFACE

BF_STUB(BfMat *, MatDenseRealGetView, BfMat *)
BF_STUB(void, MatDenseRealDelete, BfMat **)
BF_STUB(BfMat *, MatDenseRealEmptyLike, BfMat const *, BfSize, BfSize)
BF_STUB(BfMat *, MatDenseRealZerosLike, BfMat const *, BfSize, BfSize)

BfMatType bfMatDenseRealGetType(BfMat const *mat) {
  (void)mat;
  return BF_MAT_TYPE_DENSE_REAL;
}

BF_STUB(bool, MatDenseRealInstanceOf, BfMat const *, BfMatType)
BF_STUB(BfSize, MatDenseRealNumBytes, BfMat const *)
BF_STUB(void, MatDenseRealSave, BfMat const *, char const *)

void bfMatDenseRealPrint(FILE *fp, BfMat const *mat) {
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

BF_STUB(BfMat *, MatDenseRealGetRowRange, BfMat *, BfSize, BfSize)
BF_STUB(BfMat *, MatDenseRealGetColRange, BfMat *, BfSize, BfSize)
BF_STUB(void, MatDenseRealSetRowRange, BfMat *, BfSize, BfSize, BfMat const *)
BF_STUB(BfMat *, MatDenseRealRowDists, BfMat const *, BfMat const *)
BF_STUB(BfMat *, MatDenseRealColDists, BfMat const *, BfMat const *)
BF_STUB(void, MatDenseRealScaleCols, BfMat *, BfMat const *)
BF_STUB(BfMat *, MatDenseRealSumCols, BfMat const *)
BF_STUB(void, MatDenseRealAddInplace, BfMat *, BfMat const *)
BF_STUB(void, MatDenseRealAddDiag, BfMat *, BfMat const *)
BF_STUB(BfMat *, MatDenseRealMul, BfMat const *, BfMat const *)
BF_STUB(void, MatDenseRealMulInplace, BfMat *, BfMat const *)
BF_STUB(BfMat *, MatDenseRealSolve, BfMat const *, BfMat const *)
BF_STUB(BfMat *, MatDenseRealLstSq, BfMat const *, BfMat const *)

/** Upcasting: */

BfMat *bfMatDenseRealToMat(BfMatDenseReal *matDenseReal) {
  return &matDenseReal->super;
}

BfMat const *bfMatDenseRealConstToMatConst(BfMatDenseReal const *matDenseReal) {
  return &matDenseReal->super;
}

/** Downcasting: */

BfMatDenseReal const *bfMatConstToMatDenseRealConst(BfMat const *mat) {
  if (!bfMatInstanceOf(mat, BF_MAT_TYPE_DENSE_REAL)) {
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
