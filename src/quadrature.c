#include <bf/quadrature.h>

#include <assert.h>

#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mat_block_dense.h>
#include <bf/mat_coo_complex.h>
#include <bf/mat_dense_complex.h>

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
apply_KR_correction_complex(BfMatDenseComplex *mat, BfSize m, BfSize order,
                            BfKernelComplex K, BfPtr *aux) {
  BfReal const *w_KR = get_w_KR(order);
  BfComplex *rowptr;
  for (BfSize i = 0; i < m; ++i) {
    rowptr = mat->data + i*mat->rowStride;
    for (BfSize p = 0, j; p < order; ++p) {
      j = (i + p + 1) % m;
      *(rowptr + j*mat->colStride) += w_KR[p]*K(i, j, aux);
      j = ((i + m) - p - 1) % m;
      *(rowptr + j*mat->colStride) += w_KR[p]*K(i, j, aux);
    }
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

  switch (bfMatGetType(mat)) {
  case BF_TYPE_MAT_DENSE_COMPLEX:
    apply_KR_correction_complex(bfMatToMatDenseComplex(mat), m, order, K, aux);
    break;
  default:
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);
  }

  END_ERROR_HANDLING() {}
}

BfMat *bf_get_KR_corr_spmat(BfSize order, BfSize m, BfKernelComplex K, BfPtr *aux) {
  BEGIN_ERROR_HANDLING();

  BfMatCooComplex *corr = bfMatCooComplexNew();
  HANDLE_ERROR();

  bfMatCooComplexInitEmpty(corr, m, m, 2*m*order);
  HANDLE_ERROR();

  BfReal const *w_KR = get_w_KR(order);

  BfSize j, k = 0;
  for (BfSize i = 0; i < m; ++i) {
    for (BfSize p = 0; p < order; ++p) {
      j = (i + p + 1) % m;
      corr->rowInd[k] = i;
      corr->colInd[k] = j;
      corr->value[k++] = w_KR[p]*K(i, j, aux);

      j = (i - p - 1) % m;
      corr->rowInd[k] = i;
      corr->colInd[k] = j;
      corr->value[k++] = w_KR[p]*K(i, j, aux);
    }
  }

  if (k != corr->numElts)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  END_ERROR_HANDLING()
    bfMatCooComplexDeinitAndDealloc(&corr);

  return bfMatCooComplexToMat(corr);
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

  BfSize m = bfMatGetNumRows(mat);
  BfSize n = bfMatGetNumCols(mat);

  if (m != n)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (m < 2*order + 1)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  corr = bf_get_KR_corr_spmat(order, m, K, aux);
  bfMatPermuteRows(corr, perm);
  bfMatPermuteCols(corr, perm);

  bfMatAddInplace(mat, corr);

  END_ERROR_HANDLING() {}

  bfMatDelete(&corr);
}
