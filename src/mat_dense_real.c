#include <bf/mat_dense_real.h>

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <openblas/lapacke.h>

#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mat_diag_real.h>
#include <bf/vec_real.h>

/** Interface: Mat */

#define INTERFACE BF_INTERFACE_Mat
BF_DEFINE_VTABLE(Mat, MatDenseReal)
#undef INTERFACE

BfMat *bfMatDenseRealCopy(BfMat const *mat) {
  BEGIN_ERROR_HANDLING();

  BfMatDense const *matDense = NULL;
  BfMatDense *matDenseCopy = NULL;
  BfMatDenseReal const *matDenseReal = NULL;
  BfMatDenseReal *copy = NULL;

  matDense = bfMatConstToMatDenseConst(mat);
  HANDLE_ERROR();

  matDenseReal = bfMatConstToMatDenseRealConst(mat);
  HANDLE_ERROR();

  BfSize m = bfMatGetNumRows(mat);
  BfSize n = bfMatGetNumCols(mat);

  copy = bfMatDenseRealNew();
  HANDLE_ERROR();

  bfMatDenseRealInit(copy, m, n);
  HANDLE_ERROR();

  matDenseCopy = bfMatDenseRealToMatDense(copy);

  BfSize rowStride = bfMatDenseGetRowStride(matDense);
  BfSize colStride = bfMatDenseGetColStride(matDense);

  BfSize rowStrideCopy = bfMatDenseGetRowStride(matDenseCopy);
  BfSize colStrideCopy = bfMatDenseGetColStride(matDenseCopy);

  for (BfSize i = 0; i < m; ++i) {
    BfReal *outPtr = copy->data + i*rowStrideCopy;
    BfReal const *inPtr = matDenseReal->data + i*rowStride;
    for (BfSize j = 0; j < n; ++j) {
      *outPtr = *inPtr;
      outPtr += colStrideCopy;
      inPtr += colStride;
    }
  }

  END_ERROR_HANDLING()
    bfMatDenseRealDeinitAndDealloc(&copy);

  return bfMatDenseRealToMat(copy);
}

BF_STUB(BfMat *, MatDenseRealGetView, BfMat *)
BF_STUB(BfVec *, MatDenseRealGetRowCopy, BfMat const *, BfSize)
BF_STUB(BfVec *, MatDenseRealGetRowView, BfMat *, BfSize)

BfVec *bfMatDenseRealGetColView(BfMat *mat, BfSize j) {
  BEGIN_ERROR_HANDLING();

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

  END_ERROR_HANDLING()
    bfVecRealDeinitAndDealloc(&colView);

  return bfVecRealToVec(colView);
}

BF_STUB(BfVec *, MatDenseRealGetColRangeView, BfMat *, BfSize, BfSize, BfSize)

void bfMatDenseRealDelete(BfMat **mat) {
  bfMatDenseRealDeinitAndDealloc((BfMatDenseReal **)mat);
}

BF_STUB(BfMat *, MatDenseRealEmptyLike, BfMat const *, BfSize, BfSize)
BF_STUB(BfMat *, MatDenseRealZerosLike, BfMat const *, BfSize, BfSize)

BfType bfMatDenseRealGetType(BfMat const *mat) {
  (void)mat;
  return BF_TYPE_MAT_DENSE_REAL;
}

BF_STUB(bool, MatDenseRealInstanceOf, BfMat const *, BfType)
BF_STUB(BfSize, MatDenseRealNumBytes, BfMat const *)

void bfMatDenseRealSave(BfMat const *mat, char const *path) {
  BEGIN_ERROR_HANDLING();

  BfMatDenseReal const *matDenseReal = bfMatConstToMatDenseRealConst(mat);

  FILE *fp = fopen(path, "w");
  if (fp == NULL)
    RAISE_ERROR(BF_ERROR_FILE_ERROR);

  BfSize numElts = mat->numRows*mat->numCols;

  fwrite(matDenseReal->data, sizeof(BfReal), numElts, fp);
  /* TODO: error-handling */

  END_ERROR_HANDLING() {}

  fclose(fp);
}

void bfMatDenseRealPrint(BfMat const *mat, FILE *fp) {
  BEGIN_ERROR_HANDLING();

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

  END_ERROR_HANDLING() {}
}

BF_STUB(void, MatDenseRealSetCol, BfMat *, BfSize, BfVec const *)
BF_STUB(void, MatDenseRealSetColRange, BfMat *, BfSize, BfSize, BfSize, BfVec const *)
BF_STUB(BfMat *, MatDenseRealGetRowRange, BfMat *, BfSize, BfSize)
BF_STUB(BfMat *, MatDenseRealGetColRange, BfMat *, BfSize, BfSize)
BF_STUB(BfMat *, MatDenseRealGetRowRangeCopy, BfMat const *, BfSize, BfSize)
BF_STUB(BfMat *, MatDenseRealGetColRangeCopy, BfMat const *, BfSize, BfSize)
BF_STUB(void, MatDenseRealSetRowRange, BfMat *, BfSize, BfSize, BfMat const *)

void bfMatDenseRealPermuteRows(BfMat *mat, BfPerm const *perm) {
  BEGIN_ERROR_HANDLING();

  BfSize m = bfMatGetNumRows(mat);
  BfSize n = bfMatGetNumCols(mat);

  BfMatDense *matDense = NULL;
  BfMatDense *matDensePerm = NULL;
  BfMatDenseReal *matDenseReal = NULL;
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

  END_ERROR_HANDLING() {}

  bfMatDelete(&matPerm);
}

BF_STUB(void, MatDenseRealPermuteCols, BfMat *, BfPerm const *)
BF_STUB(BfVec *, MatDenseRealRowDists, BfMat const *, BfMat const *)
BF_STUB(BfVec *, MatDenseRealColDists, BfMat const *, BfMat const *)
BF_STUB(BfVec *, MatDenseRealColDots, BfMat const *, BfMat const *)
BF_STUB(BfVec *, MatDenseRealColNorms, BfMat const *)
BF_STUB(void, MatDenseRealScaleRows, BfMat *, BfVec const *)
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
BF_STUB(BfVec *, MatDenseRealForwardSolveVec, BfMat const *, BfVec const *)
BF_STUB(BfVec *, MatDenseRealBackwardSolveVec, BfMat const *, BfVec const *)
BF_STUB(bool, MatDenseRealIsZero, BfMat const *)
BF_STUB(void, MatDenseRealNegate, BfMat *)
BF_STUB(BfMat *, MatDenseRealToType, BfMat const *, BfType)
BF_STUB(BfMat *, MatDenseRealCholesky, BfMat const *)

/** Interface: MatDense */

static BfMatDenseVtable matDenseVtable = {
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

  bfMatDenseInit(&mat->super, &MatVtbl, &matDenseVtable, numRows, numCols, numCols, 1);
  HANDLE_ERROR();

  mat->data = malloc(numRows*numCols*sizeof(BfReal));
  if (mat->data == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING()
    bfMatDenseDeinit(&mat->super);
}

void bfMatDenseDeinit(BfMatDense *matDense) {
  bfMatDeinit(&matDense->super);
}

void bfMatDenseRealInitWithValue(BfMatDenseReal *matDenseReal, BfSize numRows,
                                 BfSize numCols, BfReal fillValue) {
  BEGIN_ERROR_HANDLING();

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

  END_ERROR_HANDLING()
    bfMatDenseRealDeinit(matDenseReal);
}

void bfMatDenseRealDeinit(BfMatDenseReal *matDenseReal) {
  BfMat *mat = bfMatDenseRealToMat(matDenseReal);

  if (!(mat->props & BF_MAT_PROPS_VIEW))
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

void bfMatDenseRealSvd(BfMatDenseReal const *mat, BfMatDenseReal *U,
                          BfMatDiagReal *S, BfMatDenseReal *VH) {
  BEGIN_ERROR_HANDLING();

  BfMat const *super = bfMatDenseRealConstToMatConst(mat);

  BfSize m = super->numRows;
  BfSize n = super->numCols;

  BfReal *superb = NULL;

  /* dgesvd will overwrite A, so allocate space for a copy */
  BfReal *dataCopy = malloc(m*n*sizeof(BfReal));
  if (dataCopy == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  /* copy contents of A */
  memcpy(dataCopy, mat->data, m*n*sizeof(BfReal));

  /* output array which contains information about superdiagonal
   * elements which didn't converge
   *
   * more info here: tinyurl.com/2p8f5ev3 */
  superb = malloc((((m < n) ? m : n) - 1)*sizeof(BfReal));
  if (superb == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  /* compute the SVD */
  lapack_int info = LAPACKE_dgesvd(
    LAPACK_ROW_MAJOR, 'S', 'S', m, n, dataCopy, m, S->data, U->data, m,
    VH->data, n, superb);

  /* check for invalid arguments */
  if (info < 0)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  /* check for errors */
  if (info > 0)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  END_ERROR_HANDLING() {}

  bfMatDenseRealToMat(U)->props |= BF_MAT_PROPS_ORTHO;
  bfMatDenseRealToMat(VH)->props |= BF_MAT_PROPS_ORTHO;

  free(dataCopy);
  free(superb);
}
