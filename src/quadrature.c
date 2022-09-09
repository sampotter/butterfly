#include <bf/quadrature.h>

#include <assert.h>

#include <bf/error.h>
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

static void apply_KR_correction_zmat(BfMatDenseComplex *mat, BfSize m, BfSize order) {
  BfReal const *w_KR = get_w_KR(order);
  BfComplex *rowptr;
  for (BfSize i = 0; i < m; ++i) {
    rowptr = mat->data + i*mat->rowStride;
    for (BfSize p = 0; p < order; ++p) {
      *(rowptr + ((i + p + 1) % m)*mat->colStride) *= 1 + w_KR[p];
      *(rowptr + ((i - p - 1) % m)*mat->colStride) *= 1 + w_KR[p];
    }
  }
}

void bf_apply_KR_correction(BfMat *mat, BfSize order) {
  if (order != 2 && order != 6 && order != 10) {
    bfSetError(BF_ERROR_INVALID_ARGUMENTS);
    return;
  }

  BfSize m = bfMatGetNumRows(mat);
  BfSize n = bfMatGetNumCols(mat);
  if (m != n) {
    bfSetError(BF_ERROR_INVALID_ARGUMENTS);
    return;
  }
  if (m < 2*order + 1) {
    bfSetError(BF_ERROR_INVALID_ARGUMENTS);
    return;
  }

  switch (bfMatGetType(mat)) {
  case BF_MAT_TYPE_DENSE_COMPLEX:
    apply_KR_correction_zmat(bfMatToMatDenseComplex(mat), m, order);
    break;
  default:
    bfSetError(BF_ERROR_NOT_IMPLEMENTED);
    return;
  }
}
