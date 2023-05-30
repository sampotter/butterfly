#include <bf/mat_coo_complex.h>

#include <bf/assert.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mat_dense_complex.h>
#include <bf/mem.h>
#include <bf/size_array.h>
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
  .AddInplace = (__typeof__(&bfMatAddInplace))bfMatCooComplexAddInplace,
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
  BF_ERROR_BEGIN();

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

  BF_ERROR_END()
    bfMatCooComplexDeinitAndDealloc(&rowRange);

  return bfMatCooComplexToMat(rowRange);
}

BfMat *bfMatCooComplexGetColRangeCopy(BfMat const *mat, BfSize j0, BfSize j1) {
  BF_ERROR_BEGIN();

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

  BF_ERROR_END()
    bfMatCooComplexDeinitAndDealloc(&colRange);

  return bfMatCooComplexToMat(colRange);
}

void bfMatCooComplexPermuteRows(BfMat *mat, BfPerm const *perm) {
  BF_ERROR_BEGIN();

  BfSize *rowIndPerm = NULL;

  BfMatCooComplex *matCooComplex = bfMatToMatCooComplex(mat);
  HANDLE_ERROR();

  BfSize numElts = matCooComplex->numElts;

  if (perm->size != mat->numRows)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  rowIndPerm = bfMemAlloc(numElts, sizeof(BfSize));
  if (rowIndPerm == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BfPerm revPerm = bfPermGetReversePerm(perm);
  HANDLE_ERROR();

  /* Permute the row indices */
  for (BfSize i = 0; i < numElts; ++i)
    rowIndPerm[i] = revPerm.index[matCooComplex->rowInd[i]];

  /* Free old row indices and replace with permuted ones */
  bfMemFree(matCooComplex->rowInd);
  matCooComplex->rowInd = rowIndPerm;

  BF_ERROR_END() {
    bfMemFree(rowIndPerm);
  }

  bfPermDeinit(&revPerm);
}

void bfMatCooComplexPermuteCols(BfMat *mat, BfPerm const *perm) {
  BF_ERROR_BEGIN();

  BfSize *colIndPerm = NULL;

  BfMatCooComplex *matCooComplex = bfMatToMatCooComplex(mat);
  HANDLE_ERROR();

  BfSize numElts = matCooComplex->numElts;

  if (perm->size != mat->numCols)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  colIndPerm = bfMemAlloc(numElts, sizeof(BfSize));
  if (colIndPerm == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BfPerm revPerm = bfPermGetReversePerm(perm);
  HANDLE_ERROR();

  /* Permute the column indices */
  for (BfSize i = 0; i < numElts; ++i)
    colIndPerm[i] = revPerm.index[matCooComplex->colInd[i]];

  /* Free old column indices and replace with permuted ones */
  bfMemFree(matCooComplex->colInd);
  matCooComplex->colInd = colIndPerm;

  BF_ERROR_END() {
    bfMemFree(colIndPerm);
  }

  bfPermDeinit(&revPerm);
}

static BfMat *mul_denseComplex(BfMat const *mat, BfMat const *otherMat) {
  BF_ERROR_BEGIN();

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

  BF_ERROR_END()
    bfMatDelete(&result);

  return result;
}

static void expand(BfMatCooComplex *matCooComplex) {
  BF_ERROR_BEGIN();

  BfSize newCapacity = 2*matCooComplex->capacity;

  BfSize *rowInd = NULL, *oldRowInd = NULL;
  BfSize *colInd = NULL, *oldColInd = NULL;
  BfComplex *value = NULL, *oldValue = NULL;

  rowInd = bfMemAlloc(newCapacity, sizeof(BfSize));
  if (rowInd == NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  bfMemCopy(matCooComplex->rowInd, matCooComplex->numElts, sizeof(BfSize), rowInd);
  oldRowInd = matCooComplex->rowInd;
  matCooComplex->rowInd = rowInd;

  colInd = bfMemAlloc(newCapacity, sizeof(BfSize));
  if (colInd == NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  bfMemCopy(matCooComplex->colInd, matCooComplex->numElts, sizeof(BfSize), colInd);
  oldColInd = matCooComplex->colInd;
  matCooComplex->colInd = colInd;

  value = bfMemAlloc(newCapacity, sizeof(BfComplex));
  if (value == NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  bfMemCopy(matCooComplex->value, matCooComplex->numElts, sizeof(BfComplex), value);
  oldValue = matCooComplex->value;
  matCooComplex->value = value;

  matCooComplex->capacity = newCapacity;

  BF_ERROR_END() {}

  bfMemFree(oldRowInd);
  bfMemFree(oldColInd);
  bfMemFree(oldValue);
}

static int compare(BfSize k0, BfSize k1, void *aux) {
  BfMatCooComplex *matCooComplex = (BfMatCooComplex *)aux;

  BfSize i0 = matCooComplex->rowInd[k0];
  BfSize j0 = matCooComplex->colInd[k0];

  BfSize i1 = matCooComplex->rowInd[k1];
  BfSize j1 = matCooComplex->colInd[k1];

  if (i0 < i1)
    return -1;
  else if (j0 < j1)
    return -1;
  else
    return j0 == j1 ? 0 : 1;
}

static void addInplace_cooComplex(BfMatCooComplex *matCooComplex,
                                  BfMatCooComplex const *otherMatCooComplex) {
  BF_ERROR_BEGIN();

  BfSizeArray *indexArray = bfSizeArrayNewIota(matCooComplex->numElts);
  HANDLE_ERROR();

  bfSizeArraySort(indexArray, compare, matCooComplex);

  BfSizeArray *otherIndexArray = bfSizeArrayNewIota(otherMatCooComplex->numElts);
  HANDLE_ERROR();

  bfSizeArraySort(otherIndexArray, compare, (BfMatCooComplex *)otherMatCooComplex);

  BfSize i0, i1, j0, j1;

  BfSize n0 = matCooComplex->numElts;
  BfSize n1 = otherMatCooComplex->numElts;

  BfSize k0 = 0, pk0;
  for (BfSize k1 = 0; k1 < n1; ++k1) {
    BfSize pk1 = bfSizeArrayGet(otherIndexArray, k1);

    i1 = otherMatCooComplex->rowInd[pk1];
    j1 = otherMatCooComplex->colInd[pk1];

    if (k0 < n0) {
      /* Seek k0 forward until we hit or pass k1: */
    inc:
      pk0 = bfSizeArrayGet(indexArray, k0++);
      i0 = matCooComplex->rowInd[pk0];
      j0 = matCooComplex->colInd[pk0];
      if (k0 < n0 && (i0 < i1 || (i0 == i1 && j0 < j1)))
        goto inc;

      /* If we hit k1, accumulate into k0: */
      if (i0 == i1 && j0 == j1) {
        matCooComplex->value[pk0] += otherMatCooComplex->value[pk1];
        continue;
      }
    }

    if (matCooComplex->numElts == matCooComplex->capacity) {
      expand(matCooComplex);
      HANDLE_ERROR();
    }

    matCooComplex->rowInd[matCooComplex->numElts] = i1;
    matCooComplex->colInd[matCooComplex->numElts] = j1;
    matCooComplex->value[matCooComplex->numElts] = otherMatCooComplex->value[pk1];

    ++matCooComplex->numElts;
  }

  BF_ERROR_END() {
    BF_DIE();
  }

  bfSizeArrayDeinitAndDealloc(&indexArray);
  bfSizeArrayDeinitAndDealloc(&otherIndexArray);
}

void bfMatCooComplexAddInplace(BfMatCooComplex *matCooComplex, BfMat const *otherMat) {
  switch (bfMatGetType(otherMat)) {
  case BF_TYPE_MAT_COO_COMPLEX:
    addInplace_cooComplex(matCooComplex, bfMatConstToMatCooComplexConst(otherMat));
    return;
  default:
    bfSetError(BF_ERROR_NOT_IMPLEMENTED);
  }
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
  BF_ERROR_BEGIN();

  BfMatCooComplex *mat = bfMemAlloc(1, sizeof(BfMatCooComplex));
  if (mat == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BF_ERROR_END() {}

  return mat;
}

void bfMatCooComplexInitEmpty(BfMatCooComplex *mat, BfSize numRows,
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

  mat->value = bfMemAlloc(numElts, sizeof(BfComplex));
  if (mat->value == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  mat->capacity = numElts;

  BF_ERROR_END() {
    bfMemFree(mat->rowInd);
    bfMemFree(mat->colInd);
    bfMemFree(mat->value);
    bfMatDeinit(&mat->super);
  }
}

void bfMatCooComplexDeinit(BfMatCooComplex *mat) {
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


void bfMatCooComplexDealloc(BfMatCooComplex **mat) {
  bfMemFree(*mat);
  *mat = NULL;
}

void bfMatCooComplexDeinitAndDealloc(BfMatCooComplex **mat) {
  bfMatCooComplexDeinit(*mat);
  bfMatCooComplexDealloc(mat);
}
