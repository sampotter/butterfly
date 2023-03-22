#include <bf/mat_coo_complex.h>

#include <stdlib.h>

#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mat_dense_complex.h>
#include <bf/vec_complex.h>

/** Interface: Mat */

static BfMatVtable MAT_VTABLE = {
  .Delete = (__typeof__(&bfMatCooComplexDelete))bfMatCooComplexDelete,
  .GetType = (__typeof__(&bfMatCooComplexGetType))bfMatCooComplexGetType,
  .GetNumRows = (__typeof__(&bfMatCooComplexGetNumRows))bfMatCooComplexGetNumRows,
  .GetNumCols = (__typeof__(&bfMatCooComplexGetNumCols))bfMatCooComplexGetNumCols,
  .GetRowRangeCopy = (__typeof__(&bfMatCooComplexGetRowRangeCopy))bfMatCooComplexGetRowRangeCopy,
  .GetColRangeCopy = (__typeof__(&bfMatCooComplexGetColRangeCopy))bfMatCooComplexGetColRangeCopy,
  .PermuteRows = (__typeof__(&bfMatCooComplexPermuteRows))bfMatCooComplexPermuteRows,
  .PermuteCols = (__typeof__(&bfMatCooComplexPermuteCols))bfMatCooComplexPermuteCols,
  .Mul = (__typeof__(&bfMatCooComplexMul))bfMatCooComplexMul,
  .IsZero = (__typeof__(&bfMatCooComplexIsZero))bfMatCooComplexIsZero,
};

void bfMatCooComplexDelete(BfMat **mat) {
  bfMatCooComplexDeinitAndDealloc((BfMatCooComplex **)mat);
}

BfType bfMatCooComplexGetType(BfMat const *mat) {
  (void)mat;
  return BF_TYPE_MAT_COO_COMPLEX;
}

BfSize bfMatCooComplexGetNumRows(BfMat const *mat) {
  if (bfMatGetType(mat) != BF_TYPE_MAT_COO_COMPLEX) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return BF_SIZE_BAD_VALUE;
  } else {
    return mat->numRows;
  }
}

BfSize bfMatCooComplexGetNumCols(BfMat const *mat) {
  if (bfMatGetType(mat) != BF_TYPE_MAT_COO_COMPLEX) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return BF_SIZE_BAD_VALUE;
  } else {
    return mat->numCols;
  }
}

BfMat *bfMatCooComplexGetRowRangeCopy(BfMat const *mat, BfSize i0, BfSize i1) {
  BEGIN_ERROR_HANDLING();

  if (i0 > i1)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (i1 > bfMatGetNumRows(mat))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfMatCooComplex const *matCooComplex = bfMatConstToMatCooComplexConst(mat);
  HANDLE_ERROR();

  /* Count the number of elements in the row interval [i0, i1) */
  BfSize numElts = 0;
  for (BfSize k = 0; k < matCooComplex->numElts; ++k) {
    BfSize rowInd = matCooComplex->rowInd[k];
    if (i0 <= rowInd && rowInd < i1)
      ++numElts;
  }

  BfMatCooComplex *rowRange = bfMatCooComplexNew();
  HANDLE_ERROR();

  bfMatCooComplexInitEmpty(rowRange, i1 - i0, bfMatGetNumCols(mat), numElts);
  HANDLE_ERROR();

  /* Copy over the elements in [i0, i1) */
  BfSize l = 0;
  for (BfSize k = 0; k < matCooComplex->numElts; ++k) {
    BfSize rowInd = matCooComplex->rowInd[k];
    if (i0 <= rowInd && rowInd < i1) {
      rowRange->rowInd[l] = matCooComplex->rowInd[k] - i0; /* Offset */
      rowRange->colInd[l] = matCooComplex->colInd[k];
      rowRange->value[l++] = matCooComplex->value[k];
    }
  }

  END_ERROR_HANDLING()
    bfMatCooComplexDeinitAndDealloc(&rowRange);

  return bfMatCooComplexToMat(rowRange);
}

BfMat *bfMatCooComplexGetColRangeCopy(BfMat const *mat, BfSize j0, BfSize j1) {
  BEGIN_ERROR_HANDLING();

  if (j0 > j1)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (j1 > bfMatGetNumCols(mat))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfMatCooComplex const *matCooComplex = bfMatConstToMatCooComplexConst(mat);
  HANDLE_ERROR();

  /* Count the number of elements in the col interval [j0, j1) */
  BfSize numElts = 0;
  for (BfSize k = 0; k < matCooComplex->numElts; ++k) {
    BfSize colInd = matCooComplex->colInd[k];
    if (j0 <= colInd && colInd < j1)
      ++numElts;
  }

  BfMatCooComplex *colRange = bfMatCooComplexNew();
  HANDLE_ERROR();

  bfMatCooComplexInitEmpty(colRange, bfMatGetNumRows(mat), j1 - j0, numElts);
  HANDLE_ERROR();

  /* Copy over the elements in [j0, j1) */
  BfSize l = 0;
  for (BfSize k = 0; k < matCooComplex->numElts; ++k) {
    BfSize colInd = matCooComplex->colInd[k];
    if (j0 <= colInd && colInd < j1) {
      colRange->rowInd[l] = matCooComplex->rowInd[k];
      colRange->colInd[l] = matCooComplex->colInd[k] - j0; /* Offset */
      colRange->value[l++] = matCooComplex->value[k];
    }
  }

  END_ERROR_HANDLING()
    bfMatCooComplexDeinitAndDealloc(&colRange);

  return bfMatCooComplexToMat(colRange);
}

void bfMatCooComplexPermuteRows(BfMat *mat, BfPerm const *perm) {
  BEGIN_ERROR_HANDLING();

  BfMatCooComplex *matCooComplex = bfMatToMatCooComplex(mat);
  HANDLE_ERROR();

  BfSize numElts = matCooComplex->numElts;

  if (perm->size != mat->numRows)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize *rowIndPerm = malloc(numElts*sizeof(BfSize));
  if (rowIndPerm == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BfPerm revPerm = bfPermGetReversePerm(perm);
  HANDLE_ERROR();

  /* Permute the row indices */
  for (BfSize i = 0; i < numElts; ++i)
    rowIndPerm[i] = revPerm.index[matCooComplex->rowInd[i]];

  /* Free old row indices and replace with permuted ones */
  free(matCooComplex->rowInd);
  matCooComplex->rowInd = rowIndPerm;

  END_ERROR_HANDLING() {
    free(rowIndPerm);
  }

  bfPermDeinit(&revPerm);
}

void bfMatCooComplexPermuteCols(BfMat *mat, BfPerm const *perm) {
  BEGIN_ERROR_HANDLING();

  BfMatCooComplex *matCooComplex = bfMatToMatCooComplex(mat);
  HANDLE_ERROR();

  BfSize numElts = matCooComplex->numElts;

  if (perm->size != mat->numCols)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize *colIndPerm = malloc(numElts*sizeof(BfSize));
  if (colIndPerm == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BfPerm revPerm = bfPermGetReversePerm(perm);
  HANDLE_ERROR();

  /* Permute the column indices */
  for (BfSize i = 0; i < numElts; ++i)
    colIndPerm[i] = revPerm.index[matCooComplex->colInd[i]];

  /* Free old column indices and replace with permuted ones */
  free(matCooComplex->colInd);
  matCooComplex->colInd = colIndPerm;

  END_ERROR_HANDLING() {
    free(colIndPerm);
  }

  bfPermDeinit(&revPerm);
}

static BfMat *mul_denseComplex(BfMat const *mat, BfMat const *otherMat) {
  BEGIN_ERROR_HANDLING();

  BfMatCooComplex const *matCooComplex = NULL;
  BfMatDenseComplex const *matDenseComplex = NULL;
  BfMat *result = NULL;

  BfSize m, n, p;

  matCooComplex = bfMatConstToMatCooComplexConst(mat);
  HANDLE_ERROR();

  matDenseComplex = bfMatConstToMatDenseComplexConst(otherMat);
  HANDLE_ERROR();

  m = bfMatGetNumRows(mat);
  n = bfMatGetNumCols(mat);
  p = bfMatGetNumCols(otherMat);

  if (n != bfMatGetNumRows(otherMat))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  result = bfMatZerosLike(otherMat, m, p);
  HANDLE_ERROR();

  for (BfSize k = 0; k < matCooComplex->numElts; ++k) {
    BfSize i = matCooComplex->rowInd[k];
    BfSize j = matCooComplex->colInd[k];
    BfComplex z = matCooComplex->value[k];

    BfVec *rowVec = bfMatGetRowView(result, i);
    BfVecComplex *rowVecComplex = bfVecToVecComplex(rowVec);

    BfComplex const *inPtr = matDenseComplex->data + j*matDenseComplex->rowStride;
    BfComplex *outPtr = rowVecComplex->data;

    for (BfSize l = 0; l < p; ++l) {
      *outPtr = *inPtr;
      *outPtr *= z;
      inPtr += matDenseComplex->colStride;
      outPtr += rowVecComplex->stride;
    }

    bfVecDelete(&rowVec);
  }

  END_ERROR_HANDLING()
    bfMatDelete(&result);

  return result;
}

BfMat *bfMatCooComplexMul(BfMat const *mat, BfMat const *otherMat) {
  switch (bfMatGetType(otherMat)) {
  case BF_TYPE_MAT_DENSE_COMPLEX:
    return mul_denseComplex(mat, otherMat);
  default:
    bfSetError(BF_ERROR_NOT_IMPLEMENTED);
    return NULL;
  }
}

bool bfMatCooComplexIsZero(BfMat const *mat) {
  return bfMatConstToMatCooComplexConst(mat)->numElts == 0;
}

/** Upcasting: */

BfMat *bfMatCooComplexToMat(BfMatCooComplex *matCooComplex) {
  return &matCooComplex->super;
}

BfMat const *bfMatCooComplexConstToMatConst(BfMatCooComplex const *matCooComplex) {
  return &matCooComplex->super;
}

/** Downcasting: */

BfMatCooComplex *bfMatToMatCooComplex(BfMat *mat) {
  if (!bfMatInstanceOf(mat, BF_TYPE_MAT_COO_COMPLEX)) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return NULL;
  } else {
    return (BfMatCooComplex *)mat;
  }
}

BfMatCooComplex const *bfMatConstToMatCooComplexConst(BfMat const *mat) {
  if (!bfMatInstanceOf(mat, BF_TYPE_MAT_COO_COMPLEX)) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return NULL;
  } else {
    return (BfMatCooComplex const *)mat;
  }
}

/** Implementation: MatCooComplex */

BfMatCooComplex *bfMatCooComplexNew() {
  BEGIN_ERROR_HANDLING();

  BfMatCooComplex *mat = malloc(sizeof(BfMatCooComplex));
  if (mat == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING() {}

  return mat;
}

void bfMatCooComplexInitEmpty(BfMatCooComplex *mat, BfSize numRows,
                              BfSize numCols, BfSize numElts) {
  BEGIN_ERROR_HANDLING();

  bfMatInit(&mat->super, &MAT_VTABLE, numRows, numCols);
  HANDLE_ERROR();

  mat->numElts = numElts;

  mat->rowInd = malloc(numElts*sizeof(BfSize));
  if (mat->rowInd == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  mat->colInd = malloc(numElts*sizeof(BfSize));
  if (mat->colInd == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  mat->value = malloc(numElts*sizeof(BfComplex));
  if (mat->value == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING() {
    free(mat->rowInd);
    free(mat->colInd);
    free(mat->value);
    bfMatDeinit(&mat->super);
  }
}

void bfMatCooComplexDeinit(BfMatCooComplex *mat) {
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


void bfMatCooComplexDealloc(BfMatCooComplex **mat) {
  free(*mat);
  *mat = NULL;
}

void bfMatCooComplexDeinitAndDealloc(BfMatCooComplex **mat) {
  bfMatCooComplexDeinit(*mat);
  bfMatCooComplexDealloc(mat);
}
