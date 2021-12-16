#include "mat.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <cblas.h>
#include <lapacke.h>

#include "rand.h"

BfMat bfGetUninitializedMat() {
  return (BfMat) {
    .dtype = BF_DTYPE_VOID,
    .props = BF_MAT_PROP_NONE,
    .numRows = BF_SIZE_BAD_VALUE,
    .numCols = BF_SIZE_BAD_VALUE,
    .rowStride = BF_SIZE_BAD_VALUE,
    .colStride = BF_SIZE_BAD_VALUE,
    .data = 0x0
  };
}

bool bfMatIsInitialized(BfMat const *A) {
  return A->data != NULL;
}

BfSize bfMatSize(BfMat const *A)
{
  assert(A->dtype != BF_DTYPE_MAT);

  return A->numRows*A->numCols;
}

enum BfError bfMatNumBytes(BfMat const *A, BfSize *nbytes)
{
  assert(A->dtype != BF_DTYPE_MAT);

  enum BfError error = BF_ERROR_NO_ERROR;

  BfSize dtype_size = bfDtypeSize(A->dtype);

  BfSize m = A->numRows, n = A->numCols;

  BfSize num_entries = A->props & BF_MAT_PROP_DIAGONAL ?
    (m < n ? m : n) :
    m*n;

  *nbytes = dtype_size*num_entries;

  return error;
}

BfSize bfMatNumRows(BfMat const *A) {
  return bfMatIsTransposed(A) ? A->numCols : A->numRows;
}

BfSize bfMatNumCols(BfMat const *A) {
  return bfMatIsTransposed(A) ? A->numRows : A->numCols;
}

bool bfMatIsAligned(BfMat const *A) {
  BfSize dtype_size = bfDtypeSize(A->dtype);

  return ((BfSize)A->data % dtype_size == 0)
      && (A->rowStride % dtype_size == 0)
      && (A->colStride % dtype_size == 0);
}

BfSize bfMatRowStride(BfMat const *A) {
  assert(A->dtype != BF_DTYPE_MAT);
  assert(!(A->props & BF_MAT_PROP_DIAGONAL));

  return bfMatIsTransposed(A) ? A->colStride : A->rowStride;
}

BfSize bfMatColStride(BfMat const *A) {
  assert(A->dtype != BF_DTYPE_MAT);
  assert(!(A->props & BF_MAT_PROP_DIAGONAL));

  return bfMatIsTransposed(A) ? A->rowStride : A->colStride;
}

enum BfError bfFreeMat(BfMat *A)
{
  assert(A->dtype != BF_DTYPE_MAT);

  free(A->data);

  return BF_ERROR_NO_ERROR;
}

enum BfError
bfInitEmptyMat(BfMat *A, enum BfDtypes dtype, enum BfMatProps props,
               BfSize numRows, BfSize numCols)
{
  enum BfError error = BF_ERROR_NO_ERROR;

  A->dtype = dtype;
  A->props = props;
  A->numRows = numRows;
  A->numCols = numCols;

  BfSize dtype_size = bfDtypeSize(A->dtype);

  /* when we set the row stride, if this is a diagonal matrix, we just
   * set the row stride to zero, which makes indexing simpler since we
   * just store each diagonal entry contiguously in memory */
  A->rowStride = A->props & BF_MAT_PROP_DIAGONAL ? 0 : dtype_size*numCols;

  A->colStride = dtype_size;

  A->data = malloc(dtype_size*numRows*numCols);

  return error;
}

enum BfError
bfSaveMat(BfMat const *A, char const *path)
{
  assert(A->dtype != BF_DTYPE_MAT);

  enum BfError error = BF_ERROR_NO_ERROR;

  FILE *fp = fopen(path, "w");
  if (fp == NULL) {
    error = BF_ERROR_RUNTIME_ERROR;
    goto cleanup;
  }

  BfSize dtype_size = bfDtypeSize(A->dtype);
  BfSize m = A->numRows, n = A->numCols;
  BfSize num_entries = A->props & BF_MAT_PROP_DIAGONAL ?
    (m < n ? m : n) :
    m*n;

  fwrite(A->data, dtype_size, num_entries, fp);

cleanup:
  fclose(fp);
  return error;
}

enum BfError
bfFillMatRandn(BfMat *A)
{
  assert(A->dtype != BF_DTYPE_MAT);

  switch (A->dtype) {
  case BF_DTYPE_REAL:
    bfRandn(bfMatSize(A), A->data);
    return BF_ERROR_NO_ERROR;
  case BF_DTYPE_COMPLEX:
    bfRandn(2*bfMatSize(A), A->data);
    return BF_ERROR_NO_ERROR;
  default:
    if (bfDtypeIsValid(A->dtype))
      return BF_ERROR_NOT_IMPLEMENTED;
    else
      return BF_ERROR_RUNTIME_ERROR;
  }
}

static bool validIndex(BfMat const *A, BfSize i, BfSize j) {
  return i < A->numRows && j < A->numCols;
}

static BfSize offset(BfMat const *A, BfSize i, BfSize j) {
  if (A->props & BF_MAT_PROP_DIAGONAL)
    assert(i == j);

  return i*A->rowStride + j*A->colStride;
}

static void
getDiagonalMatElt(BfMat const *A, BfSize i, BfSize j, BfPtr ptr) {
  assert(A->dtype != BF_DTYPE_MAT);

  if (i == j) {
    BfSize size = bfDtypeSize(A->dtype);
    memcpy(ptr, (BfByte *)A->data + i*size, size);
  } else {
    // TODO: we should make sure to return a zero matrix with a
    // conformable shape here
    bfGetDtypeZero(A->dtype, ptr);
  }
}

enum BfError
bfGetMatElt(BfMat const *A, BfSize i, BfSize j, BfPtr ptr) {
  assert(A->dtype != BF_DTYPE_MAT);

  if (!validIndex(A, i, j))
    return BF_ERROR_OUT_OF_RANGE;

  if (A->props & BF_MAT_PROP_DIAGONAL) {
    getDiagonalMatElt(A, i, j, ptr);
  } else {
    BfSize dtype_size = bfDtypeSize(A->dtype);
    memcpy(ptr, (BfByte *)A->data + offset(A, i, j), dtype_size);
  }

  return BF_ERROR_NO_ERROR;
}

enum BfError
getDiagonalMatEltPtr(BfMat const *A, BfSize i, BfSize j, BfPtr *ptr) {
  if (i != j)
    return BF_ERROR_INVALID_ARGUMENTS;

  BfSize size = bfDtypeSize(A->dtype);

  *ptr = (BfByte *)A->data + i*size;

  return BF_ERROR_NO_ERROR;
}

enum BfError
bfGetMatEltPtr(BfMat const *A, BfSize i, BfSize j, BfPtr *ptr) {
  if (!validIndex(A, i, j))
    return BF_ERROR_OUT_OF_RANGE;

  if (A->props & BF_MAT_PROP_DIAGONAL)
    return getDiagonalMatEltPtr(A, i, j, ptr);

  *ptr = (BfByte *)A->data + offset(A, i, j);

  return BF_ERROR_NO_ERROR;
}

BfVec bfGetMatRow(BfMat const *A, BfSize i) {
  assert(A->dtype != BF_DTYPE_MAT);
  assert(!(A->props & BF_MAT_PROP_DIAGONAL));

  return (BfVec) {
    .dtype = A->dtype,
    .props = BF_VEC_PROP_VIEW,
    .size = bfMatNumCols(A),
    .stride = bfMatColStride(A),
    .data = (BfByte *)A->data + i*bfMatRowStride(A)
  };
}

enum BfError
bfGetRowPtr(BfMat const *A, BfSize i, BfPtr *ptr)
{
  assert(A->dtype != BF_DTYPE_MAT);
  assert(!(A->props & BF_MAT_PROP_DIAGONAL));

  if (i >= bfMatNumRows(A))
    return BF_ERROR_OUT_OF_RANGE;

  *ptr = (BfByte *)A->data + i*bfMatRowStride(A);

  return BF_ERROR_NO_ERROR;
}

enum BfError
bfSetMatRow(BfMat *A, BfSize i, void const *data)
{
  assert(A->dtype != BF_DTYPE_MAT);
  assert(!(A->props & BF_MAT_PROP_DIAGONAL));

  BfSize m = A->numRows, n = A->numCols;

  if (i >= m)
    return BF_ERROR_INVALID_ARGUMENTS;

  BfSize nbytes = bfDtypeSize(A->dtype);

  BfSize row_nbytes = nbytes*n;

  memcpy((BfByte *)A->data + row_nbytes*i, data, row_nbytes);

  return BF_ERROR_NO_ERROR;
}

enum BfError
bfCopyMatRow(BfMat const *A, BfSize i, BfMat *B, BfSize j)
{
  assert(A->dtype != BF_DTYPE_MAT);
  assert(!(A->props & BF_MAT_PROP_DIAGONAL));

  BfSize m = bfMatNumCols(A);

  if (m != bfMatNumCols(B))
    return BF_ERROR_BAD_SHAPE;

  if (i >= bfMatNumRows(A))
    return BF_ERROR_OUT_OF_RANGE;

  if (j >= bfMatNumRows(B))
    return BF_ERROR_OUT_OF_RANGE;

  if (A->dtype != B->dtype)
    return BF_ERROR_INVALID_ARGUMENTS;

  BfByte *A_row_ptr, *B_row_ptr;
  bfGetRowPtr(A, i, (BfPtr *)&A_row_ptr);
  bfGetRowPtr(B, j, (BfPtr *)&B_row_ptr);

  BfSize A_col_stride = bfMatColStride(A);
  BfSize B_col_stride = bfMatColStride(B);

  // TODO: this sucks
  if (A->dtype == BF_DTYPE_COMPLEX) {
    bool shouldConj = (A->props & BF_MAT_PROP_CONJ_TRANS)
                    ^ (B->props & BF_MAT_PROP_CONJ_TRANS);

    for (BfSize k = 0; k < m; ++k) {
      *(BfComplex *)B_row_ptr = shouldConj ?
        conj(*(BfComplex *)A_row_ptr) : *(BfComplex *)A_row_ptr;

      A_row_ptr += A_col_stride;
      B_row_ptr += B_col_stride;
    }
  } else {
    assert(false); // TODO: not implemented
  }

  return BF_ERROR_NO_ERROR;
}

BfMat bfGetMatRowRange(BfMat const *A, BfSize i0, BfSize i1) {
  assert(A->dtype != BF_DTYPE_MAT);
  assert(!(A->props & BF_MAT_PROP_DIAGONAL));

  assert(!bfMatIsTransposed(A));

  BfMat A_rows = *A;

  A_rows.props |= BF_MAT_PROP_VIEW;
  if (A_rows.props & BF_MAT_PROP_UNITARY) {
    A_rows.props ^= BF_MAT_PROP_UNITARY;
    A_rows.props ^= BF_MAT_PROP_SEMI_UNITARY;
  }

  A_rows.numRows = i1 - i0;
  A_rows.data = (BfByte *)A_rows.data + A->rowStride*i0;

  return A_rows;
}

BfMat bfGetMatColRange(BfMat const *A, BfSize j0, BfSize j1) {
  assert(A->dtype != BF_DTYPE_MAT);
  assert(!(A->props & BF_MAT_PROP_DIAGONAL));

  assert(!bfMatIsTransposed(A));

  BfMat A_cols = *A;

  A_cols.props |= BF_MAT_PROP_VIEW;

  if (A_cols.props & BF_MAT_PROP_UNITARY) {
    A_cols.props ^= BF_MAT_PROP_UNITARY;
    A_cols.props ^= BF_MAT_PROP_SEMI_UNITARY;
  }

  A_cols.numCols = j1 - j0;
  A_cols.data = (BfByte *)A_cols.data + A->colStride*j0;

  return A_cols;
}

BfMat bfGetMatContSubblock(BfMat const *A, BfSize i0, BfSize i1,
                           BfSize j0, BfSize j1) {
  assert(A->dtype != BF_DTYPE_MAT);

  // TODO: if i0 != j0, we need to return a block matrix with zero
  // subblocks except for a single diagonal subblock... this is
  // tedious so just disallow this for now
  if (A->dtype & BF_MAT_PROP_DIAGONAL)
    assert(i0 == j0);

  BfMat A_subblock = *A;

  A_subblock.props |= BF_MAT_PROP_VIEW;

  A_subblock.numRows = i1 - i0;
  A_subblock.numCols = j1 - j0;

  A_subblock.data = (BfByte *)A_subblock.data + offset(A, i0, j0);

  return A_subblock;
}

bool bfMatIsTransposed(BfMat const *A) {
  return A->props & (BF_MAT_PROP_TRANS | BF_MAT_PROP_CONJ_TRANS);
}

BfMat bfConjTrans(BfMat const *A) {
  assert(A->dtype != BF_DTYPE_MAT);
  assert(!(A->props & BF_MAT_PROP_DIAGONAL));

  BfMat AH = *A;

  /* make AH a view of A */
  AH.props |= BF_MAT_PROP_VIEW;

  /* make AH the Hermitian transpose of A */
  AH.props ^= BF_MAT_PROP_CONJ_TRANS;

  return AH;
}

static enum CBLAS_TRANSPOSE getCblasTranspose(BfMat const *A) {
  assert(A->dtype != BF_DTYPE_MAT);
  assert(!(A->props & BF_MAT_PROP_DIAGONAL));

  if (A->props & BF_MAT_PROP_CONJ_TRANS)
    return CblasConjTrans;
  else if (A->props & BF_MAT_PROP_TRANS)
    return CblasTrans;
  else
    return CblasNoTrans;
}

static BfSize getLeadingDimension(BfMat const *A) {
  assert(A->dtype != BF_DTYPE_MAT);
  assert(!(A->props & BF_MAT_PROP_DIAGONAL));

  BfSize dtype_size = bfDtypeSize(A->dtype);
  return bfMatIsTransposed(A) ?
    bfMatColStride(A)/dtype_size :
    bfMatRowStride(A)/dtype_size;
}

static void
complexComplexMatMul(BfMat const *A, BfMat const *B, BfMat *C)
{
  assert(A->dtype == BF_DTYPE_COMPLEX);
  assert(B->dtype == BF_DTYPE_COMPLEX);

  assert(!(A->props & BF_MAT_PROP_DIAGONAL));
  assert(!(B->props & BF_MAT_PROP_DIAGONAL));

  assert(!bfMatIsInitialized(C));

  BfComplex alpha = 1, beta = 0;

  enum CBLAS_TRANSPOSE A_trans = getCblasTranspose(A);
  enum CBLAS_TRANSPOSE B_trans = getCblasTranspose(B);

  BfSize m = bfMatNumRows(A);
  BfSize n = bfMatNumCols(B);
  BfSize k = bfMatNumCols(A);
  assert(k == bfMatNumRows(B));
  bfInitEmptyMat(C, BF_DTYPE_COMPLEX, BF_MAT_PROP_NONE, m, n);

  BfSize lda = getLeadingDimension(A);
  BfSize ldb = getLeadingDimension(B);
  BfSize ldc = getLeadingDimension(C);

  cblas_zgemm(CblasRowMajor, A_trans, B_trans, m, n, k,
              &alpha, A->data, lda, B->data, ldb, &beta, C->data, ldc);
}

enum BfError
bfMatMul(BfMat const *A, BfMat const *B, BfMat *C)
{
  assert(A->dtype != BF_DTYPE_MAT);
  assert(B->dtype != BF_DTYPE_MAT);

  assert(!(A->props & BF_MAT_PROP_DIAGONAL));

  assert(bfMatIsAligned(A));
  assert(bfMatIsAligned(B));

  if (bfMatNumCols(A) != bfMatNumRows(B))
    return BF_ERROR_INVALID_ARGUMENTS;

  if (A->dtype == BF_DTYPE_COMPLEX && B->dtype == BF_DTYPE_COMPLEX)
    complexComplexMatMul(A, B, C);
  else
    return BF_ERROR_NOT_IMPLEMENTED;

  return BF_ERROR_NO_ERROR;
}

enum BfError
realComplexMatSolve(BfMat const *A, BfMat const *B, BfMat *C)
{
  assert(A->dtype == BF_DTYPE_REAL);
  assert(B->dtype == BF_DTYPE_COMPLEX);

  BfSize m = bfMatNumCols(A);
  BfSize k = bfMatNumRows(B);
  BfSize n = bfMatNumCols(B);

  bfInitEmptyMat(C, BF_DTYPE_COMPLEX, BF_MAT_PROP_NONE, m, n);

  if (A->props & BF_MAT_PROP_DIAGONAL) {
    BfReal scale;
    BfVec row;
    for (BfSize i = 0; i < k; ++i) {
      bfCopyMatRow(B, i, C, i);
      row = bfGetMatRow(C, i);
      bfGetMatElt(A, i, i, &scale);
      bfVecScaleByReal(&row, 1/scale);
    }
  }

  return BF_ERROR_NO_ERROR;
}

enum BfError
complexComplexMatSolve(BfMat const *A, BfMat const *B, BfMat *C)
{
  assert(!(A->props & BF_MAT_PROP_DIAGONAL));
  assert(!(A->props & BF_MAT_PROP_SEMI_UNITARY));

  if (A->props & BF_MAT_PROP_UNITARY) {
    BfMat AH = *A;
    AH.props ^= (AH.props & BF_MAT_PROP_CONJ_TRANS);
    complexComplexMatMul(&AH, B, C);
    return BF_ERROR_NO_ERROR;
  }

  return BF_ERROR_NOT_IMPLEMENTED;
}

/* Solve the linear system B = A*C for C. Analogous to A\B = C in
 * MATLAB. */
enum BfError
bfMatSolve(BfMat const *A, BfMat const *B, BfMat *C)
{
  assert(A->dtype != BF_DTYPE_MAT);
  assert(B->dtype != BF_DTYPE_MAT);

  if (bfMatNumRows(A) != bfMatNumRows(B))
    return BF_ERROR_INVALID_ARGUMENTS;

  if (A->dtype == BF_DTYPE_REAL && B->dtype == BF_DTYPE_COMPLEX)
    return realComplexMatSolve(A, B, C);
  else if (A->dtype == BF_DTYPE_COMPLEX && B->dtype == BF_DTYPE_COMPLEX)
    return complexComplexMatSolve(A, B, C);
  else
    return BF_ERROR_NOT_IMPLEMENTED;
}

enum BfError
initEmptySvdMatsComplex(BfMat const *A, BfMat *U, BfMat *S, BfMat *Vt)
{
  assert(A->dtype != BF_DTYPE_MAT);
  assert(!(A->props & BF_MAT_PROP_DIAGONAL));

  enum BfError error = BF_ERROR_NO_ERROR;

  BfSize m = A->numRows;
  BfSize n = A->numCols;
  BfSize p = m < n ? m : n;

  /* init U */
  error = bfInitEmptyMat(U, BF_DTYPE_COMPLEX, BF_MAT_PROP_NONE, m, p);
  if (error)
    return error;

  /* init S */
  error = bfInitEmptyMat(S, BF_DTYPE_REAL, BF_MAT_PROP_DIAGONAL, p, p);
  if (error) {
    bfFreeMat(U);
    return error;
  }

  /* init Vt */
  error = bfInitEmptyMat(Vt, BF_DTYPE_COMPLEX, BF_MAT_PROP_NONE, p, n);
  if (error) {
    bfFreeMat(U);
    bfFreeMat(S);
    return error;
  }

  return error;
}

enum BfError
bfInitEmptySvdMats(BfMat const *A, BfMat *U, BfMat *S, BfMat *Vt)
{
  assert(A->dtype != BF_DTYPE_MAT);
  assert(!(A->props & BF_MAT_PROP_DIAGONAL));

  switch (A->dtype) {
  case BF_DTYPE_COMPLEX:
    return initEmptySvdMatsComplex(A, U, S, Vt);
  default:
    return BF_ERROR_NOT_IMPLEMENTED;
  }
}

enum BfError
computeMatSvdComplex(BfMat const *A, BfMat *U, BfMat *S, BfMat *Vt)
{
  assert(A->dtype != BF_DTYPE_MAT);
  assert(!(A->props & BF_MAT_PROP_DIAGONAL));

  enum BfError error = BF_ERROR_NO_ERROR;

  lapack_int m = A->numRows;
  lapack_int n = A->numCols;

  /* zgesvd will overwrite A, so make a copy of it here for input */
  BfComplex *A_data_copy = malloc(m*n*sizeof(BfComplex));
  memcpy(A_data_copy, A->data, m*n*sizeof(BfComplex));

  /* output array which contains information about superdiagonal
   * elements which didn't converge
   *
   * more info here: tinyurl.com/2p8f5ev3 */
  BfReal *superb = malloc((((m < n) ? m : n) - 1)*sizeof(BfReal));

  /* compute the SVD */
  lapack_int info = LAPACKE_zgesvd(
    LAPACK_ROW_MAJOR, 'S', 'S', m, n, A_data_copy, m, S->data, U->data, m,
    Vt->data, n, superb);

  /* check for invalid arguments */
  if (info < 0) {
    error = BF_ERROR_INVALID_ARGUMENTS;
    goto cleanup;
  }

  /* check for errors */
  if (info > 0) {
    error = BF_ERROR_RUNTIME_ERROR;
    goto cleanup;
  }

  U->props |= BF_MAT_PROP_UNITARY;
  S->props |= BF_MAT_PROP_DIAGONAL;
  Vt->props |= BF_MAT_PROP_UNITARY;

cleanup:
  free(A_data_copy);
  free(superb);

  return error;
}

enum BfError
bfComputeMatSvd(BfMat const *A, BfMat *U, BfMat *S, BfMat *Vt)
{
  assert(A->dtype != BF_DTYPE_MAT);
  assert(!(A->props & BF_MAT_PROP_DIAGONAL));

  bfInitEmptySvdMats(A, U, S, Vt);

  switch (A->dtype) {
  case BF_DTYPE_COMPLEX:
    return computeMatSvdComplex(A, U, S, Vt);
  default:
    return BF_ERROR_NOT_IMPLEMENTED;
  }
}

enum BfError
bfComputePinv(BfMat const *A, BfReal atol, BfReal rtol, BfMat *pinv)
{
  assert(A->dtype != BF_DTYPE_MAT);
  assert(!(A->props & BF_MAT_PROP_DIAGONAL));

  enum BfError error;

  /* compute SVD of A */

  BfMat U = bfGetUninitializedMat();
  BfMat S = bfGetUninitializedMat();
  BfMat VH = bfGetUninitializedMat();
  error = bfComputeMatSvd(A, &U, &S, &VH);
  if (error)
    goto cleanup;

  /* compute tolerance and compute number of terms to
   * retain in pseudoinverse */

  BfSize m = A->numRows, n = A->numCols;
  BfSize k_max = m < n ? m : n;

  BfReal *sigma = S.data;

  BfReal tol = rtol*sigma[0] + atol;

  BfSize k;
  for (k = 0; k < k_max; ++k)
    if (sigma[k] < tol)
      break;

  /* get subblocks of truncated SVD */

  BfMat UkH = bfGetMatColRange(&U, 0, k);
  UkH = bfConjTrans(&UkH);

  BfMat Sk = bfGetMatContSubblock(&S, 0, k, 0, k);

  BfMat Vk = bfGetMatRowRange(&VH, 0, k);
  Vk = bfConjTrans(&Vk);

  /* compute pseudonverse from truncated SVD */

  BfMat WkH = bfGetUninitializedMat();
  error = bfMatSolve(&Sk, &UkH, &WkH);
  if (error)
    goto cleanup;

  error = bfMatMul(&Vk, &WkH, pinv);

cleanup:
  bfFreeMat(&U);
  bfFreeMat(&S);
  bfFreeMat(&VH);
  bfFreeMat(&WkH);

  return error;
}

enum BfError
bfMatLstSq(BfMat const *A, BfMat const *B, BfMat *C)
{
  assert(!(A->props & BF_MAT_PROP_DIAGONAL));

  enum BfError error = BF_ERROR_NO_ERROR;

  BfReal const atol = BF_EPS_MACH;

  BfSize m = A->numRows;
  BfSize n = A->numCols;
  BfReal const rtol = (m > n ? m : n)*BF_EPS_MACH;

  BfMat A_pinv = bfGetUninitializedMat();
  error = bfComputePinv(A, atol, rtol, &A_pinv);
  if (error)
    goto cleanup;

  bfSaveMat(&A_pinv, "Z_eq_pinv.bin");

  error = bfMatMul(&A_pinv, B, C);
  if (error)
    goto cleanup;

cleanup:
  bfFreeMat(&A_pinv);

  return error;
}
