#include <bf/mat_csr_real.h>

#include <stdlib.h>
#include <string.h>

#include <bf/error.h>
#include <bf/error_macros.h>

/** Interface: Mat */

#define INTERFACE BF_INTERFACE_Mat
BF_DEFINE_VTABLE(Mat, MatCsrReal)
#undef INTERFACE

BF_STUB(BfMat *, MatCsrRealCopy, BfMat const *)
BF_STUB(BfMat *, MatCsrRealGetView, BfMat *)
BF_STUB(BfVec *, MatCsrRealGetRowCopy, BfMat const *, BfSize)
BF_STUB(BfVec *, MatCsrRealGetRowView, BfMat *, BfSize)
BF_STUB(BfVec *, MatCsrRealGetColView, BfMat *, BfSize)
BF_STUB(BfVec *, MatCsrRealGetColRangeView, BfMat *, BfSize, BfSize, BfSize)

void bfMatCsrRealDelete(BfMat **mat) {
  bfMatCsrRealDeinitAndDealloc((BfMatCsrReal **)mat);
}

BF_STUB(BfMat *, MatCsrRealEmptyLike, BfMat const *, BfSize, BfSize)
BF_STUB(BfMat *, MatCsrRealZerosLike, BfMat const *, BfSize, BfSize)

BfType bfMatCsrRealGetType(BfMat const *mat) {
  (void)mat;
  return BF_TYPE_MAT_CSR_REAL;
}

BF_STUB(BfSize, MatCsrRealNumBytes, BfMat const *)
BF_STUB(void, MatCsrRealSave, BfMat const *, char const *)
BF_STUB(void, MatCsrRealPrint, BfMat const *, FILE *)

BfSize bfMatCsrRealGetNumRows(BfMat const *mat) {
  if (bfMatGetType(mat) != BF_TYPE_MAT_CSR_REAL) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return BF_SIZE_BAD_VALUE;
  } else {
    return mat->numRows;
  }
}

BfSize bfMatCsrRealGetNumCols(BfMat const *mat) {
  if (bfMatGetType(mat) != BF_TYPE_MAT_CSR_REAL) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return BF_SIZE_BAD_VALUE;
  } else {
    return mat->numCols;
  }
}

BF_STUB(void, MatCsrRealSetRow, BfMat *, BfSize, BfVec const *)
BF_STUB(void, MatCsrRealSetCol, BfMat *, BfSize, BfVec const *)
BF_STUB(void, MatCsrRealSetColRange, BfMat *, BfSize, BfSize, BfSize, BfVec const *)
BF_STUB(BfMat *, MatCsrRealGetRowRange, BfMat *, BfSize, BfSize)
BF_STUB(BfMat *, MatCsrRealGetColRange, BfMat *, BfSize, BfSize)
BF_STUB(BfMat *, MatCsrRealGetRowRangeCopy, BfMat const *, BfSize, BfSize)
BF_STUB(BfMat *, MatCsrRealGetColRangeCopy, BfMat const *, BfSize, BfSize)
BF_STUB(void, MatCsrRealSetRowRange, BfMat *, BfSize, BfSize, BfMat const *)
BF_STUB(void, MatCsrRealPermuteRows, BfMat *, BfPerm const *)
BF_STUB(void, MatCsrRealPermuteCols, BfMat *, BfPerm const *)
BF_STUB(BfVec *, MatCsrRealRowDists, BfMat const *, BfMat const *)
BF_STUB(BfVec *, MatCsrRealColDists, BfMat const *, BfMat const *)
BF_STUB(BfVec *, MatCsrRealColDots, BfMat const *, BfMat const *)
BF_STUB(BfVec *, MatCsrRealColNorms, BfMat const *)
BF_STUB(void, MatCsrRealScaleRows, BfMat *, BfVec const *)
BF_STUB(void, MatCsrRealScaleCols, BfMat *, BfVec const *)
BF_STUB(BfVec *, MatCsrRealSumCols, BfMat const *)
BF_STUB(void, MatCsrRealAddInplace, BfMat *, BfMat const *)
BF_STUB(void, MatCsrRealAddDiag, BfMat *, BfMat const *)
BF_STUB(BfMat *, MatCsrRealSub, BfMat const *, BfMat const *)
BF_STUB(void, MatCsrRealSubInplace, BfMat *, BfMat const *)
BF_STUB(BfMat *, MatCsrRealMul, BfMat const *, BfMat const *)
BF_STUB(BfVec *, MatCsrRealMulVec, BfMat const *, BfVec const *)
BF_STUB(void, MatCsrRealMulInplace, BfMat *, BfMat const *)
BF_STUB(BfMat *, MatCsrRealSolveLU, BfMat const *, BfMat const *)
BF_STUB(BfMat *, MatCsrRealLstSq, BfMat const *, BfMat const *)
BF_STUB(bool, MatCsrRealIsUpperTri, BfMat const *)
BF_STUB(BfVec *, MatCsrRealForwardSolveVec, BfMat const *, BfVec const *)
BF_STUB(BfVec *, MatCsrRealBackwardSolveVec, BfMat const *, BfVec const *)

bool bfMatCsrRealIsZero(BfMat const *mat) {
  return bfMatConstToMatCsrRealConst(mat)->rowptr[mat->numRows] == 0;
}

BF_STUB(void, MatCsrRealNegate, BfMat *)
BF_STUB(BfMat *, MatCsrRealToType, BfMat const *, BfType)

/** Upcasting: */

BfMat *bfMatCsrRealToMat(BfMatCsrReal *matCsrReal) {
  return &matCsrReal->super;
}

BfMat const *bfMatCsrRealConstToMatConst(BfMatCsrReal const *matCsrReal) {
  return &matCsrReal->super;
}

/** Downcasting: */

BfMatCsrReal *bfMatToMatCsrReal(BfMat *mat) {
  if (!bfMatInstanceOf(mat, BF_TYPE_MAT_CSR_REAL)) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return NULL;
  } else {
    return (BfMatCsrReal *)mat;
  }
}

BfMatCsrReal const *bfMatConstToMatCsrRealConst(BfMat const *mat) {
  if (!bfMatInstanceOf(mat, BF_TYPE_MAT_CSR_REAL)) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return NULL;
  } else {
    return (BfMatCsrReal const *)mat;
  }
}

/** Implementation: MatCsrReal */

BfMatCsrReal *bfMatCsrRealNew() {
  BEGIN_ERROR_HANDLING();

  BfMatCsrReal *mat = malloc(sizeof(BfMatCsrReal));
  if (mat == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING() {}

  return mat;
}

void bfMatCsrRealInit(BfMatCsrReal *mat, BfSize numRows, BfSize numCols,
                      BfSize const *rowptr, BfSize const *colind,
                      BfReal const *data) {
  BEGIN_ERROR_HANDLING();

  BfSize nnz = rowptr[numRows];

  bfMatInit(&mat->super, &MatVtbl, numRows, numCols);
  HANDLE_ERROR();

  mat->rowptr = malloc((numRows + 1)*sizeof(BfSize));
  if (mat->rowptr == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  mat->colind = malloc(nnz*sizeof(BfSize));
  if (mat->colind == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  mat->data = malloc(nnz*sizeof(BfReal));
  if (mat->data == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  memcpy(mat->rowptr, rowptr, (numRows + 1)*sizeof(BfSize));
  memcpy(mat->colind, colind, nnz*sizeof(BfSize));
  memcpy(mat->data, data, nnz*sizeof(BfSize));

  END_ERROR_HANDLING()
    bfMatCsrRealDeinit(mat);
}

void bfMatCsrRealDeinit(BfMatCsrReal *mat) {
  if (!(mat->super.props & BF_MAT_PROPS_VIEW)) {
    free(mat->rowptr);
    free(mat->colind);
    free(mat->data);
  }

  mat->rowptr = NULL;
  mat->colind = NULL;
  mat->data = NULL;

  bfMatDeinit(&mat->super);
}

void bfMatCsrRealDealloc(BfMatCsrReal **mat) {
  free(*mat);
  *mat = NULL;
}

void bfMatCsrRealDeinitAndDealloc(BfMatCsrReal **mat) {
  bfMatCsrRealDeinit(*mat);
  bfMatCsrRealDealloc(mat);
}

void bfMatCsrRealDump(BfMatCsrReal const *matCsrReal, char const *rowptrPath,
                      char const *colindPath, char const *dataPath) {
  BEGIN_ERROR_HANDLING();

  BfSize numRows = matCsrReal->super.numRows;
  BfSize nnz = matCsrReal->rowptr[numRows];

  FILE *fp;

  /* Save `rowptr` to disk: */

  fp = fopen(rowptrPath, "w");
  if (fp == NULL)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);
  fwrite(matCsrReal->rowptr, sizeof(BfSize), numRows + 1, fp);
  fclose(fp);

  /* Save `colind` to disk: */

  fp = fopen(colindPath, "w");
  if (fp == NULL)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);
  fwrite(matCsrReal->colind, sizeof(BfSize), nnz, fp);
  fclose(fp);

  /* Save `data` to disk: */

  fp = fopen(dataPath, "w");
  if (fp == NULL)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);
  fwrite(matCsrReal->data, sizeof(BfReal), nnz, fp);
  fclose(fp);

  END_ERROR_HANDLING() {
    // TODO: delete all files
    fclose(fp);
  }
}
