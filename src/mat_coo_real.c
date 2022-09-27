#include <bf/mat_coo_real.h>

#include <stdlib.h>

#include <bf/error.h>
#include <bf/error_macros.h>

/** Interface: Mat */

#define INTERFACE BF_INTERFACE_Mat
BF_DEFINE_VTABLE(Mat, MatCooReal)
#undef INTERFACE

BF_STUB(BfMat *, MatCooRealCopy, BfMat const *)
BF_STUB(BfMat *, MatCooRealGetView, BfMat *)
BF_STUB(BfVec *, MatCooRealGetRowView, BfMat *, BfSize)
BF_STUB(BfVec *, MatCooRealGetColView, BfMat *, BfSize)
BF_STUB(BfVec *, MatCooRealGetColRangeView, BfMat *, BfSize, BfSize, BfSize)

void bfMatCooRealDelete(BfMat **mat) {
  bfMatCooRealDeinitAndDealloc((BfMatCooReal **)mat);
}

BF_STUB(BfMat *, MatCooRealEmptyLike, BfMat const *, BfSize, BfSize)
BF_STUB(BfMat *, MatCooRealZerosLike, BfMat const *, BfSize, BfSize)

BfType bfMatCooRealGetType(BfMat const *mat) {
  (void)mat;
  return BF_TYPE_MAT_COO_REAL;
}

BF_STUB(BfSize, MatCooRealNumBytes, BfMat const *)
BF_STUB(void, MatCooRealSave, BfMat const *, char const *)
BF_STUB(void, MatCooRealPrint, BfMat const *, FILE *)

BfSize bfMatCooRealGetNumRows(BfMat const *mat) {
  if (bfMatGetType(mat) != BF_TYPE_MAT_COO_REAL) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return BF_SIZE_BAD_VALUE;
  } else {
    return mat->numRows;
  }
}

BfSize bfMatCooRealGetNumCols(BfMat const *mat) {
  if (bfMatGetType(mat) != BF_TYPE_MAT_COO_REAL) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return BF_SIZE_BAD_VALUE;
  } else {
    return mat->numCols;
  }
}

BF_STUB(void, MatCooRealSetRow, BfMat *, BfSize, BfVec const *)
BF_STUB(void, MatCooRealSetCol, BfMat *, BfSize, BfVec const *)
BF_STUB(void, MatCooRealSetColRange, BfMat *, BfSize, BfSize, BfSize, BfVec const *)
BF_STUB(BfMat *, MatCooRealGetRowRange, BfMat *, BfSize, BfSize)
BF_STUB(BfMat *, MatCooRealGetColRange, BfMat *, BfSize, BfSize)

BfMat *bfMatCooRealGetRowRangeCopy(BfMat const *mat, BfSize i0, BfSize i1) {
  BEGIN_ERROR_HANDLING();

  if (i0 > i1)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (i1 > bfMatGetNumRows(mat))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfMatCooReal const *matCooReal = bfMatConstToMatCooRealConst(mat);
  HANDLE_ERROR();

  /* Count the number of elements in the row interval [i0, i1) */
  BfSize numElts = 0;
  for (BfSize k = 0; k < matCooReal->numElts; ++k) {
    BfSize rowInd = matCooReal->rowInd[k];
    if (i0 <= rowInd && rowInd < i1)
      ++numElts;
  }

  BfMatCooReal *rowRange = bfMatCooRealNew();
  HANDLE_ERROR();

  bfMatCooRealInitEmpty(rowRange, i1 - i0, bfMatGetNumCols(mat), numElts);
  HANDLE_ERROR();

  /* Copy over the elements in [i0, i1) */
  BfSize l = 0;
  for (BfSize k = 0; k < matCooReal->numElts; ++k) {
    BfSize rowInd = matCooReal->rowInd[k];
    if (i0 <= rowInd && rowInd < i1) {
      rowRange->rowInd[l] = matCooReal->rowInd[k] - i0; /* Offset */
      rowRange->colInd[l] = matCooReal->colInd[k];
      rowRange->value[l++] = matCooReal->value[k];
    }
  }

  END_ERROR_HANDLING()
    bfMatCooRealDeinitAndDealloc(&rowRange);

  return bfMatCooRealToMat(rowRange);
}

BfMat *bfMatCooRealGetColRangeCopy(BfMat const *mat, BfSize j0, BfSize j1) {
  BEGIN_ERROR_HANDLING();

  if (j0 > j1)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (j1 > bfMatGetNumCols(mat))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfMatCooReal const *matCooReal = bfMatConstToMatCooRealConst(mat);
  HANDLE_ERROR();

  /* Count the number of elements in the col interval [j0, j1) */
  BfSize numElts = 0;
  for (BfSize k = 0; k < matCooReal->numElts; ++k) {
    BfSize colInd = matCooReal->colInd[k];
    if (j0 <= colInd && colInd < j1)
      ++numElts;
  }

  BfMatCooReal *colRange = bfMatCooRealNew();
  HANDLE_ERROR();

  bfMatCooRealInitEmpty(colRange, bfMatGetNumRows(mat), j1 - j0, numElts);
  HANDLE_ERROR();

  /* Copy over the elements in [j0, j1) */
  BfSize l = 0;
  for (BfSize k = 0; k < matCooReal->numElts; ++k) {
    BfSize colInd = matCooReal->colInd[k];
    if (j0 <= colInd && colInd < j1) {
      colRange->rowInd[l] = matCooReal->rowInd[k];
      colRange->colInd[l] = matCooReal->colInd[k] - j0; /* Offset */
      colRange->value[l++] = matCooReal->value[k];
    }
  }

  END_ERROR_HANDLING()
    bfMatCooRealDeinitAndDealloc(&colRange);

  return bfMatCooRealToMat(colRange);
}

BF_STUB(void, MatCooRealSetRowRange, BfMat *, BfSize, BfSize, BfMat const *)

void bfMatCooRealPermuteRows(BfMat *mat, BfPerm const *perm) {
  BEGIN_ERROR_HANDLING();

  BfMatCooReal *matCooReal = bfMatToMatCooReal(mat);
  HANDLE_ERROR();

  BfSize numElts = matCooReal->numElts;

  if (perm->size != mat->numRows)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize *rowIndPerm = malloc(numElts*sizeof(BfSize));
  if (rowIndPerm == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BfPerm revPerm = bfPermGetReversePerm(perm);
  HANDLE_ERROR();

  /* Permute the row indices */
  for (BfSize i = 0; i < numElts; ++i)
    rowIndPerm[i] = revPerm.index[matCooReal->rowInd[i]];

  /* Free old row indices and replace with permuted ones */
  free(matCooReal->rowInd);
  matCooReal->rowInd = rowIndPerm;

  END_ERROR_HANDLING() {
    free(rowIndPerm);
  }

  bfPermDeinit(&revPerm);
}

void bfMatCooRealPermuteCols(BfMat *mat, BfPerm const *perm) {
  BEGIN_ERROR_HANDLING();

  BfMatCooReal *matCooReal = bfMatToMatCooReal(mat);
  HANDLE_ERROR();

  BfSize numElts = matCooReal->numElts;

  if (perm->size != mat->numCols)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize *colIndPerm = malloc(numElts*sizeof(BfSize));
  if (colIndPerm == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BfPerm revPerm = bfPermGetReversePerm(perm);
  HANDLE_ERROR();

  /* Permute the column indices */
  for (BfSize i = 0; i < numElts; ++i)
    colIndPerm[i] = revPerm.index[matCooReal->colInd[i]];

  /* Free old column indices and replace with permuted ones */
  free(matCooReal->colInd);
  matCooReal->colInd = colIndPerm;

  END_ERROR_HANDLING() {
    free(colIndPerm);
  }

  bfPermDeinit(&revPerm);
}

BF_STUB(BfVec *, MatCooRealRowDists, BfMat const *, BfMat const *)
BF_STUB(BfVec *, MatCooRealColDists, BfMat const *, BfMat const *)
BF_STUB(BfVec *, MatCooRealColDots, BfMat const *, BfMat const *)
BF_STUB(BfVec *, MatCooRealColNorms, BfMat const *)
BF_STUB(void, MatCooRealScaleCols, BfMat *, BfVec const *)
BF_STUB(BfVec *, MatCooRealSumCols, BfMat const *)
BF_STUB(void, MatCooRealAddInplace, BfMat *, BfMat const *)
BF_STUB(void, MatCooRealAddDiag, BfMat *, BfMat const *)
BF_STUB(BfMat *, MatCooRealSub, BfMat const *, BfMat const *)
BF_STUB(void, MatCooRealSubInplace, BfMat *, BfMat const *)
BF_STUB(BfMat *, MatCooRealMul, BfMat const *, BfMat const *)
BF_STUB(BfVec *, MatCooRealMulVec, BfMat const *, BfVec const *)
BF_STUB(void, MatCooRealMulInplace, BfMat *, BfMat const *)
BF_STUB(BfMat *, MatCooRealSolveLU, BfMat const *, BfMat const *)
BF_STUB(BfMat *, MatCooRealLstSq, BfMat const *, BfMat const *)
BF_STUB(bool, MatCooRealIsUpperTri, BfMat const *)
BF_STUB(BfVec *, MatCooRealBackwardSolveVec, BfMat const *, BfVec const *)

/** Upcasting: */

BfMat *bfMatCooRealToMat(BfMatCooReal *matCooReal) {
  return &matCooReal->super;
}

BfMat const *bfMatCooRealConstToMatConst(BfMatCooReal const *matCooReal) {
  return &matCooReal->super;
}

/** Downcasting: */

BfMatCooReal *bfMatToMatCooReal(BfMat *mat) {
  if (!bfMatInstanceOf(mat, BF_TYPE_MAT_COO_REAL)) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return NULL;
  } else {
    return (BfMatCooReal *)mat;
  }
}

BfMatCooReal const *bfMatConstToMatCooRealConst(BfMat const *mat) {
  if (!bfMatInstanceOf(mat, BF_TYPE_MAT_COO_REAL)) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return NULL;
  } else {
    return (BfMatCooReal const *)mat;
  }
}

/** Implementation: MatCooReal */

BfMatCooReal *bfMatCooRealNew() {
  BEGIN_ERROR_HANDLING();

  BfMatCooReal *mat = malloc(sizeof(BfMatCooReal));
  if (mat == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING() {}

  return mat;
}

void bfMatCooRealInitEmpty(BfMatCooReal *mat, BfSize numRows,
                              BfSize numCols, BfSize numElts) {
  BEGIN_ERROR_HANDLING();

  bfMatInit(&mat->super, &MatVtbl, numRows, numCols);
  HANDLE_ERROR();

  mat->numElts = numElts;

  mat->rowInd = malloc(numElts*sizeof(BfSize));
  if (mat->rowInd == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  mat->colInd = malloc(numElts*sizeof(BfSize));
  if (mat->colInd == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  mat->value = malloc(numElts*sizeof(BfReal));
  if (mat->value == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING() {
    free(mat->rowInd);
    free(mat->colInd);
    free(mat->value);
    bfMatDeinit(&mat->super);
  }
}

void bfMatCooRealDeinit(BfMatCooReal *mat) {
  if (!(mat->super.props & BF_MAT_PROPS_VIEW)) {
    free(mat->rowInd);
    free(mat->colInd);
    free(mat->value);
  }

  mat->numElts = BF_SIZE_BAD_VALUE;
  mat->rowInd = NULL;
  mat->colInd = NULL;
  mat->value = NULL;
}


void bfMatCooRealDealloc(BfMatCooReal **mat) {
  free(*mat);
  *mat = NULL;
}

void bfMatCooRealDeinitAndDealloc(BfMatCooReal **mat) {
  bfMatCooRealDeinit(*mat);
  bfMatCooRealDealloc(mat);
}
