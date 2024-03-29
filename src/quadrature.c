#include <bf/quadrature.h>

#include <bf/assert.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mat_block_dense.h>
#include <bf/mat_block_diag.h>
#include <bf/mat_coo_complex.h>
#include <bf/mat_dense_complex.h>
#include <bf/size_array.h>

/* Weights for 2nd order Kapur-Rokhlin quadrature. */
static BfReal const w_KR_order2[2] = {
   1.825748064736159,
  -1.325748064736159
};

/* Weights for 6th order Kapur-Rokhlin quadrature. */
static BfReal const w_KR_order6[6] = {
    4.967362978287758,
  -16.20501504859126,
   25.85153761832639,
  -22.22599466791883,
    9.930104998037539,
   -1.817995878141594
};

/* Weights for 10th order Kapur-Rokhlin quadrature. */
static BfReal const w_KR_order10[10] = {
   7.832432020568779,
  -4.565161670374749,
   1.452168846354677,
  -2.901348302886379,
   3.870862162579900,
  -3.523821383570681,
   2.172421547519342,
  -8.707796087382991,
   2.053584266072635,
  -2.166984103403823
};

static BfReal const *get_w_KR(BfSize order) {
  switch (order) {
  case 2: return w_KR_order2;
  case 6: return w_KR_order6;
  case 10: return w_KR_order10;
  default: BF_DIE();
  }
}

void bfQuadKrAccumCorrection(BfSize order, BfKernelComplex K, BfPtr *aux,
                             BfSize n, BfComplex const *x, BfReal const *h, BfComplex *y) {
  BF_ERROR_BEGIN();

  if (order != 2 && order != 6 && order != 10)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (n < 2*order + 1)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfReal const *w_KR = get_w_KR(order);

  for (BfSize i = 0; i < n; ++i) {
    for (BfSize p = 0, j; p < order; ++p) {
      j = (i + p + 1) % n;         y[i] += h[j]*w_KR[p]*K(i, j, aux)*x[j];
      j = ((i + n) - p - 1) % n;   y[i] += h[j]*w_KR[p]*K(i, j, aux)*x[j];
    }
  }

  BF_ERROR_END() {
    BF_DIE();
  }
}

static void
applyBlockCorrection_complex(BfMatDenseComplex *mat, BfSize i0, BfSize i1, BfSize order,
                             BfKernelComplex K, BfPtr *aux) {
  BfReal const *w_KR = get_w_KR(order);
  BfComplex *rowptr;
  BfSize m = i1 - i0;
  for (BfSize i = i0; i < i1; ++i) {
    rowptr = bfMatDenseComplexGetRowPtr(mat, i);
    for (BfSize p = 0, j; p < order; ++p) {
      j = ((i + p + 1 - i0) % m) + i0;
      *(rowptr + j*mat->colStride) += w_KR[p]*K(i, j, aux);
      j = (((i + m) - p - 1 - i0) % m) + i0;
      *(rowptr + j*mat->colStride) += w_KR[p]*K(i, j, aux);
    }
  }
}

static void applyBlockCorrection(BfMat *mat, BfSize i0, BfSize i1, BfSize order,
                                 BfKernelComplex K, BfPtr *aux) {
  switch (bfMatGetType(mat)) {
  case BF_TYPE_MAT_DENSE_COMPLEX:
    applyBlockCorrection_complex(bfMatToMatDenseComplex(mat), i0, i1, order, K, aux);
    break;
  default:
    bfSetError(BF_ERROR_NOT_IMPLEMENTED);
  }
}

void bfQuadKrApplyCorrection(BfMat *mat, BfSize order, BfKernelComplex K, BfPtr *aux) {
  BF_ERROR_BEGIN();

  if (order != 2 && order != 6 && order != 10)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize m = bfMatGetNumRows(mat);
  BfSize n = bfMatGetNumCols(mat);

  if (m != n)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (m < 2*order + 1)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  applyBlockCorrection(mat, 0, m, order, K, aux);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }
}

static BfMat *bf_get_KR_corr_block_spmat(BfSize order, BfSize n, BfSize i0, BfSize i1, BfKernelComplex K, BfPtr *aux) {
  BF_ERROR_BEGIN();

  if (i0 > i1)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (i1 > n)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfMatCooComplex *corr = bfMatCooComplexNew();
  HANDLE_ERROR();

  BfSize m = i1 - i0;

  bfMatCooComplexInitEmpty(corr, n, n, 2*m*order);
  HANDLE_ERROR();

  BfReal const *w_KR = get_w_KR(order);

  BfSize j, k = 0;
  for (BfSize i = i0; i < i1; ++i) {
    for (BfSize p = 0; p < order; ++p) {
      j = ((i + p + 1 - i0) % m) + i0;
      corr->rowInd[k] = i;
      corr->colInd[k] = j;
      corr->value[k++] = w_KR[p]*K(i, j, aux);

      j = (((i + m) - p - 1 - i0) % m) + i0;
      corr->rowInd[k] = i;
      corr->colInd[k] = j;
      corr->value[k++] = w_KR[p]*K(i, j, aux);
    }
  }

  if (k != 2*m*order)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  BF_ERROR_END()
    bfMatCooComplexDeinitAndDealloc(&corr);

  return bfMatCooComplexToMat(corr);
}

static BfMat *bf_get_KR_corr_spmat(BfSize order, BfSize n, BfKernelComplex K, BfPtr *aux) {
  return bf_get_KR_corr_block_spmat(order, n, 0, n, K, aux);
}

void bfQuadKrApplyCorrectionTree(BfMat *mat, BfSize order, BfTree const *tree, BfKernelComplex K, BfPtr *aux) {
  BF_ERROR_BEGIN();

  BfMat *corr = NULL;

  BfPerm const *perm = bfTreeGetPermConst(tree);

  if (order != 2 && order != 6 && order != 10)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize n = bfMatGetNumRows(mat);
  if (n != bfMatGetNumCols(mat))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (n < 2*order + 1)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  corr = bf_get_KR_corr_spmat(order, n, K, aux);
  bfMatPermuteRows(corr, perm);
  bfMatPermuteCols(corr, perm);

  bfMatAddInplace(mat, corr);

  BF_ERROR_END() {}

  bfMatDelete(&corr);
}

void bfQuadKrApplyBlockCorrection(BfMat *mat, BfSizeArray const *offsets, BfSize order, BfKernelComplex K, BfPtr *aux) {
  BF_ERROR_BEGIN();

  if (bfSizeArrayGetSize(offsets) < 2)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (!bfSizeArrayIsSorted(offsets))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize i0 = bfSizeArrayGetFirst(offsets);
  for (BfSize i = 1; i < bfSizeArrayGetSize(offsets); ++i) {
    BfSize i1 = bfSizeArrayGet(offsets, i);
    applyBlockCorrection(mat, i0, i1, order, K, aux);
    i0 = i1;
  }

  BF_ERROR_END() {
    BF_DIE();
  }
}

void bfQuadKrApplyBlockCorrectionTree(BfMat *mat, BfSizeArray const *offsets, BfSize order, BfTree const *tree, BfKernelComplex K, BfPtr *aux) {
  BF_ERROR_BEGIN();

  if (bfSizeArrayGetSize(offsets) < 2)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (!bfSizeArrayIsSorted(offsets))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize n = bfMatGetNumRows(mat);
  if (n != bfMatGetNumCols(mat))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  /** Get KR correction sparse matrix for each diagonal block
   ** corresponding to a pair of offsets: */

  BfMat *corr = NULL;

  BfSize i0 = bfSizeArrayGetFirst(offsets);
  for (BfSize i = 1; i < bfSizeArrayGetSize(offsets); ++i) {
    BfSize i1 = bfSizeArrayGet(offsets, i);
    BfMat *corrBlock = bf_get_KR_corr_block_spmat(order, n, i0, i1, K, aux);
    if (corr == NULL) {
      corr = corrBlock;
    } else {
      bfMatAddInplace(corr, corrBlock);
      bfMatDelete(&corrBlock);
    }
    i0 = i1;
  }

  /** Permute block diagonal correction matrix's rows and columns
   ** using tree: */

  BfPerm const *perm = bfTreeGetPermConst(tree);
  bfMatPermuteRows(corr, perm);
  bfMatPermuteCols(corr, perm);

  /** Apply the correction to `mat` by adding in place: */

  bfMatAddInplace(mat, corr);

  BF_ERROR_END() {
    BF_DIE();
  }

  bfMatDelete(&corr);
}
