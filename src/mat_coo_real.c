#include <bf/mat_coo_real.h>

#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mem.h>

/** Interface: Mat */

static BfMatVtable MAT_VTABLE = {
  .Delete = (__typeof__(&bfMatCooRealDelete))bfMatCooRealDelete,
  .GetType = (__typeof__(&bfMatCooRealGetType))bfMatCooRealGetType,
  .GetNumRows = (__typeof__(&bfMatCooRealGetNumRows))bfMatCooRealGetNumRows,
  .GetNumCols = (__typeof__(&bfMatCooRealGetNumCols))bfMatCooRealGetNumCols,
  .GetRowRangeCopy = (__typeof__(&bfMatCooRealGetRowRangeCopy))bfMatCooRealGetRowRangeCopy,
  .GetColRangeCopy = (__typeof__(&bfMatCooRealGetColRangeCopy))bfMatCooRealGetColRangeCopy,
  .PermuteRows = (__typeof__(&bfMatCooRealPermuteRows))bfMatCooRealPermuteRows,
  .PermuteCols = (__typeof__(&bfMatCooRealPermuteCols))bfMatCooRealPermuteCols,
  .IsZero = (__typeof__(&bfMatCooRealIsZero))bfMatCooRealIsZero,
};

void bfMatCooRealDelete(BfMat **mat) {
  bfMatCooRealDeinitAndDealloc((BfMatCooReal **)mat);
}

BfType bfMatCooRealGetType(BfMat const *mat) {
  (void)mat;
  return BF_TYPE_MAT_COO_REAL;
}

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

BfMat *bfMatCooRealGetRowRangeCopy(BfMat const *mat, BfSize i0, BfSize i1) {
  BF_ERROR_BEGIN();

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
  BF_ERROR_BEGIN();

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

void bfMatCooRealPermuteRows(BfMat *mat, BfPerm const *perm) {
  BF_ERROR_BEGIN();

  BfSize *rowIndPerm = NULL;

  BfMatCooReal *matCooReal = bfMatToMatCooReal(mat);
  HANDLE_ERROR();

  BfSize numElts = matCooReal->numElts;

  if (perm->size != mat->numRows)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  rowIndPerm = bfMemAlloc(numElts, sizeof(BfSize));
  if (rowIndPerm == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BfPerm revPerm = bfPermGetReversePerm(perm);
  HANDLE_ERROR();

  /* Permute the row indices */
  for (BfSize i = 0; i < numElts; ++i)
    rowIndPerm[i] = revPerm.index[matCooReal->rowInd[i]];

  /* Free old row indices and replace with permuted ones */
  bfMemFree(matCooReal->rowInd);
  matCooReal->rowInd = rowIndPerm;

  END_ERROR_HANDLING() {
    bfMemFree(rowIndPerm);
  }

  bfPermDeinit(&revPerm);
}

void bfMatCooRealPermuteCols(BfMat *mat, BfPerm const *perm) {
  BF_ERROR_BEGIN();

  BfSize *colIndPerm = NULL;

  BfMatCooReal *matCooReal = bfMatToMatCooReal(mat);
  HANDLE_ERROR();

  BfSize numElts = matCooReal->numElts;

  if (perm->size != mat->numCols)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  colIndPerm = bfMemAlloc(numElts, sizeof(BfSize));
  if (colIndPerm == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BfPerm revPerm = bfPermGetReversePerm(perm);
  HANDLE_ERROR();

  /* Permute the column indices */
  for (BfSize i = 0; i < numElts; ++i)
    colIndPerm[i] = revPerm.index[matCooReal->colInd[i]];

  /* Free old column indices and replace with permuted ones */
  bfMemFree(matCooReal->colInd);
  matCooReal->colInd = colIndPerm;

  END_ERROR_HANDLING() {
    bfMemFree(colIndPerm);
  }

  bfPermDeinit(&revPerm);
}

bool bfMatCooRealIsZero(BfMat const *mat) {
  return bfMatConstToMatCooRealConst(mat)->numElts == 0;
}

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
  BF_ERROR_BEGIN();

  BfMatCooReal *mat = bfMemAlloc(1, sizeof(BfMatCooReal));
  if (mat == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING() {}

  return mat;
}

void bfMatCooRealInitEmpty(BfMatCooReal *mat, BfSize numRows,
                              BfSize numCols, BfSize numElts) {
  BF_ERROR_BEGIN();

  bfMatInit(&mat->super, &MAT_VTABLE, numRows, numCols);
  HANDLE_ERROR();

  mat->numElts = numElts;

  mat->rowInd = bfMemAlloc(numElts, sizeof(BfSize));
  if (mat->rowInd == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  mat->colInd = bfMemAlloc(numElts, sizeof(BfSize));
  if (mat->colInd == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  mat->value = bfMemAlloc(numElts, sizeof(BfReal));
  if (mat->value == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING() {
    bfMemFree(mat->rowInd);
    bfMemFree(mat->colInd);
    bfMemFree(mat->value);
    bfMatDeinit(&mat->super);
  }
}

void bfMatCooRealDeinit(BfMatCooReal *mat) {
  if (!(mat->super.props & BF_MAT_PROPS_VIEW)) {
    bfMemFree(mat->rowInd);
    bfMemFree(mat->colInd);
    bfMemFree(mat->value);
  }

  mat->numElts = BF_SIZE_BAD_VALUE;
  mat->rowInd = NULL;
  mat->colInd = NULL;
  mat->value = NULL;
}

void bfMatCooRealDealloc(BfMatCooReal **mat) {
  bfMemFree(*mat);
  *mat = NULL;
}

void bfMatCooRealDeinitAndDealloc(BfMatCooReal **mat) {
  bfMatCooRealDeinit(*mat);
  bfMatCooRealDealloc(mat);
}
