#include <bf/quadrature.h>

#include <assert.h>

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
  default: assert(false);
  }
}

static void
apply_KR_correction_to_block_complex(BfMatDenseComplex *mat, BfSize i0, BfSize i1, BfSize order,
                                     BfKernelComplex K, BfPtr *aux) {
  BfReal const *w_KR = get_w_KR(order);
  BfComplex *rowptr;
  BfSize m = i1 - i0;
  for (BfSize i = i0; i < i1; ++i) {
    rowptr = mat->data + i*mat->rowStride;
    for (BfSize p = 0, j; p < order; ++p) {
      j = ((i + p + 1 - i0) % m) + i0;
      *(rowptr + j*mat->colStride) += w_KR[p]*K(i, j, aux);
      j = (((i + m) - p - 1 - i0) % m) + i0;
      *(rowptr + j*mat->colStride) += w_KR[p]*K(i, j, aux);
    }
  }
}

static void apply_KR_correction_to_block(BfMat *mat, BfSize i0, BfSize i1, BfSize order,
                                         BfKernelComplex K, BfPtr *aux) {
  switch (bfMatGetType(mat)) {
  case BF_TYPE_MAT_DENSE_COMPLEX:
    apply_KR_correction_to_block_complex(bfMatToMatDenseComplex(mat), i0, i1, order, K, aux);
    break;
  default:
    bfSetError(BF_ERROR_NOT_IMPLEMENTED);
  }
}

void bf_apply_KR_correction(BfMat *mat, BfSize order,
                            BfKernelComplex K, BfPtr *aux) {
  BEGIN_ERROR_HANDLING();

  if (order != 2 && order != 6 && order != 10)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize m = bfMatGetNumRows(mat);
  BfSize n = bfMatGetNumCols(mat);

  if (m != n)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (m < 2*order + 1)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  apply_KR_correction_to_block(mat, 0, m, order, K, aux);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {
    assert(false);
  }
}

BfMat *bf_get_KR_corr_block_spmat(BfSize order, BfSize n, BfSize i0, BfSize i1, BfKernelComplex K, BfPtr *aux) {
  BEGIN_ERROR_HANDLING();

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

  END_ERROR_HANDLING()
    bfMatCooComplexDeinitAndDealloc(&corr);

  return bfMatCooComplexToMat(corr);
}

BfMat *bf_get_KR_corr_spmat(BfSize order, BfSize n, BfKernelComplex K, BfPtr *aux) {
  return bf_get_KR_corr_block_spmat(order, n, 0, n, K, aux);
}

void bf_apply_KR_correction_quadtree(BfMat *mat, BfSize order,
                                     BfTree const *tree,
                                     BfKernelComplex K, BfPtr *aux)
{
  BEGIN_ERROR_HANDLING();

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

  END_ERROR_HANDLING() {}

  bfMatDelete(&corr);
}

void bf_apply_block_KR_correction(BfMat *mat, BfSizeArray const *offsets, BfSize order, BfKernelComplex K, BfPtr *aux) {
  BEGIN_ERROR_HANDLING();

  if (bfSizeArrayGetSize(offsets) < 2)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (!bfSizeArrayIsSorted(offsets))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize i0 = bfSizeArrayGetFirst(offsets);
  for (BfSize i = 1; i < bfSizeArrayGetSize(offsets); ++i) {
    BfSize i1 = bfSizeArrayGet(offsets, i);
    apply_KR_correction_to_block(mat, i0, i1, order, K, aux);
    i0 = i1;
  }

  END_ERROR_HANDLING() {
    assert(false);
  }
}

void bf_apply_block_KR_correction_quadtree(BfMat *mat, BfSizeArray const *offsets, BfSize order, BfTree const *tree, BfKernelComplex K, BfPtr *aux) {
  BEGIN_ERROR_HANDLING();

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

  END_ERROR_HANDLING() {
    assert(false);
  }

  bfMatDelete(&corr);
}
