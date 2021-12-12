#include "mat.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <cblas.h>
#include <lapacke.h>

#include "rand.h"

BfSize bfMatSize(BfMat const *A)
{
  assert(A->dtype != BF_DTYPE_MAT);

  return A->shape[0]*A->shape[1];
}

enum BfError bfMatNumBytes(BfMat const *A, BfSize *nbytes)
{
  assert(A->dtype != BF_DTYPE_MAT);

  enum BfError error = BF_ERROR_NO_ERROR;

  BfSize dtype_size;
  error |= bfSizeOfDtype(A->dtype, &dtype_size);
  if (error)
    return error;

  BfSize m = A->shape[0], n = A->shape[1];

  BfSize num_entries = A->props & BF_MAT_PROP_DIAGONAL ?
    (m < n ? m : n) :
    m*n;

  *nbytes = dtype_size*num_entries;

  return error;
}

BfSize bfMatNumRows(BfMat const *A) {
  assert(A->dtype != BF_DTYPE_MAT);

  return bfMatIsTransposed(A) ? A->shape[1] : A->shape[0];
}

BfSize bfMatNumCols(BfMat const *A) {
  assert(A->dtype != BF_DTYPE_MAT);

  return bfMatIsTransposed(A) ? A->shape[0] : A->shape[1];
}

bool bfMatIsAligned(BfMat const *A) {
  assert(A->dtype != BF_DTYPE_MAT);

  BfSize dtype_size = bfDtypeSize(A->dtype);
  bool ptr_aligned = (BfSize)A->data % dtype_size == 0;
  bool stride_0_aligned = A->stride[0] % dtype_size == 0;
  bool stride_1_aligned = A->stride[1] % dtype_size == 0;
  return ptr_aligned && stride_0_aligned && stride_1_aligned;
}

BfSize bfMatRowStride(BfMat const *A) {
  assert(A->dtype != BF_DTYPE_MAT);

  return bfMatIsTransposed(A) ? A->stride[1] : A->stride[0];
}

BfSize bfMatColStride(BfMat const *A) {
  return bfMatIsTransposed(A) ? A->stride[0] : A->stride[1];
}

enum BfError bfFreeMat(BfMat *A)
{
  assert(A->dtype != BF_DTYPE_MAT);

  free(A->data);

  return BF_ERROR_NO_ERROR;
}

enum BfError
bfInitEmptyMat(BfMat *A, enum BfDtypes dtype, enum BfMatProps props,
               BfSize const shape[2])
{
  assert(A->dtype != BF_DTYPE_MAT);

  enum BfError error = BF_ERROR_NO_ERROR;

  A->dtype = dtype;
  A->props = props;
  A->shape[0] = shape[0];
  A->shape[1] = shape[1];

  BfSize dtype_size;
  error |= bfSizeOfDtype(A->dtype, &dtype_size);
  if (error)
    return error;

  A->stride[0] = dtype_size*shape[1];
  A->stride[1] = dtype_size;

  A->data = malloc(dtype_size*shape[0]*shape[1]);

  return error;
}

enum BfError
bfSaveMat(BfMat const *A, char const *path)
{
  assert(A->dtype != BF_DTYPE_MAT);

  enum BfError error = BF_ERROR_NO_ERROR;

  FILE *fp = fopen(path, "w");

  BfSize dtype_size;
  error = bfSizeOfDtype(A->dtype, &dtype_size);
  if (error)
    goto cleanup;

  BfSize m = A->shape[0], n = A->shape[1];

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
  return i < A->shape[0] && j < A->shape[1];
}

static BfSize offset(BfMat const *A, BfSize i, BfSize j) {
  assert(A->props & ~BF_MAT_PROP_DIAGONAL);

  return i*A->stride[0] + j*A->stride[1];
}

static void
getDiagonalMatElt(BfMat const *A, BfSize i, BfSize j, BfPtr ptr) {
  assert(A->dtype != BF_DTYPE_MAT);

  if (i == j) {
    BfSize size = bfDtypeSize(A->dtype);
    memcpy(ptr, A->data + i*size, size);
  } else {
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
    memcpy(ptr, A->data + offset(A, i, j), dtype_size);
  }

  return BF_ERROR_NO_ERROR;
}

enum BfError
getDiagonalMatEltPtr(BfMat const *A, BfSize i, BfSize j, BfPtr *ptr) {
  if (i != j)
    return BF_ERROR_INVALID_ARGUMENTS;

  BfSize size = bfDtypeSize(A->dtype);

  *ptr = A->data + i*size;

  return BF_ERROR_NO_ERROR;
}

enum BfError
bfGetMatEltPtr(BfMat const *A, BfSize i, BfSize j, BfPtr *ptr) {
  if (!validIndex(A, i, j))
    return BF_ERROR_OUT_OF_RANGE;

  if (A->props & BF_MAT_PROP_DIAGONAL)
    return getDiagonalMatEltPtr(A, i, j, ptr);

  *ptr = A->data + offset(A, i, j);

  return BF_ERROR_NO_ERROR;
}

BfVec bfGetMatRow(BfMat const *A, BfSize i) {
  assert(A->dtype != BF_DTYPE_MAT);

  return (BfVec) {
    .dtype = A->dtype,
    .props = BF_VEC_PROP_VIEW,
    .size = bfMatNumCols(A),
    .stride = bfMatColStride(A),
    .data = A->data + i*bfMatRowStride(A)
  };
}

enum BfError
bfGetRowPtr(BfMat const *A, BfSize i, BfPtr *ptr)
{
  assert(A->dtype != BF_DTYPE_MAT);

  if (i >= bfMatNumRows(A))
    return BF_ERROR_OUT_OF_RANGE;

  *ptr = A->data + i*bfMatRowStride(A);

  return BF_ERROR_NO_ERROR;
}

enum BfError
bfSetMatRow(BfMat *A, BfSize i, void const *data)
{
  assert(A->dtype != BF_DTYPE_MAT);

  BfSize m = A->shape[0], n = A->shape[1];

  if (i >= m)
    return BF_ERROR_INVALID_ARGUMENTS;

  BfSize nbytes;
  enum BfError error = bfSizeOfDtype(A->dtype, &nbytes);
  if (error)
    return error;

  BfSize row_nbytes = nbytes*n;

  memcpy(A->data + row_nbytes*i, data, row_nbytes);

  return BF_ERROR_NO_ERROR;
}

enum BfError
bfCopyMatRow(BfMat const *A, BfSize i, BfMat *B, BfSize j)
{
  assert(A->dtype != BF_DTYPE_MAT);

  BfSize m = bfMatNumCols(A);

  if (m != bfMatNumCols(B))
    return BF_ERROR_BAD_SHAPE;

  if (i >= bfMatNumRows(A))
    return BF_ERROR_OUT_OF_RANGE;

  if (j >= bfMatNumRows(B))
    return BF_ERROR_OUT_OF_RANGE;

  BfByte *A_row_ptr, *B_row_ptr;
  bfGetRowPtr(A, i, (BfPtr *)&A_row_ptr);
  bfGetRowPtr(B, j, (BfPtr *)&B_row_ptr);

  BfSize A_col_stride = bfMatColStride(A);
  BfSize B_col_stride = bfMatColStride(B);

  for (BfSize k = 0; k < m; ++k) {
    B_row_ptr[k] = A_row_ptr[k];
    A_row_ptr += A_col_stride;
    B_row_ptr += B_col_stride;
  }

  return BF_ERROR_NO_ERROR;
}

BfMat bfGetMatRowRange(BfMat const *A, BfSize i0, BfSize i1) {
  assert(A->dtype != BF_DTYPE_MAT);

  BfMat A_rows = *A;

  A_rows.props |= BF_MAT_PROP_VIEW;
  if (A_rows.props & BF_MAT_PROP_UNITARY) {
    A_rows.props ^= BF_MAT_PROP_UNITARY;
    A_rows.props ^= BF_MAT_PROP_SEMI_UNITARY;
  }

  A_rows.shape[0] = i1 - i0;
  A_rows.data += i0*A->shape[1];

  return A_rows;
}

BfMat bfGetMatColRange(BfMat const *A, BfSize j0, BfSize j1) {
  assert(A->dtype != BF_DTYPE_MAT);

  BfMat A_cols = *A;

  A_cols.props |= BF_MAT_PROP_VIEW;

  if (A_cols.props & BF_MAT_PROP_UNITARY) {
    A_cols.props ^= BF_MAT_PROP_UNITARY;
    A_cols.props ^= BF_MAT_PROP_SEMI_UNITARY;
  }

  A_cols.shape[1] = j1 - j0;
  A_cols.data += j0;

  return A_cols;
}

BfMat bfGetMatContSubblock(BfMat const *A, BfSize i0, BfSize i1,
                           BfSize j0, BfSize j1) {
  assert(A->dtype != BF_DTYPE_MAT);

  BfMat A_subblock = *A;
  A_subblock.props |= BF_MAT_PROP_VIEW;
  A_subblock.shape[0] = i1 - i0;
  A_subblock.shape[1] = j1 - j0;
  A_subblock.data += i0*A->shape[1] + j0;
  return A_subblock;
}

bool bfMatIsTransposed(BfMat const *A) {
  assert(A->dtype != BF_DTYPE_MAT);

  return A->props & (BF_MAT_PROP_TRANS | BF_MAT_PROP_CONJ_TRANS);
}

BfMat bfConjTrans(BfMat const *A) {
  assert(A->dtype != BF_DTYPE_MAT);

  BfMat AH = *A;

  /* make AH a view of A */
  AH.props |= BF_MAT_PROP_VIEW;

  /* make AH the Hermitian transpose of A */
  AH.props ^= BF_MAT_PROP_CONJ_TRANS;

  return AH;
}

static enum CBLAS_TRANSPOSE getCblasTranspose(BfMat const *A) {
  assert(A->dtype != BF_DTYPE_MAT);

  if (A->props & BF_MAT_PROP_CONJ_TRANS)
    return CblasConjTrans;
  else if (A->props & BF_MAT_PROP_TRANS)
    return CblasTrans;
  else
    return CblasNoTrans;
}

static BfSize getLeadingDimension(BfMat const *A) {
  assert(A->dtype != BF_DTYPE_MAT);

  BfSize dtype_size = bfDtypeSize(A->dtype);
  return bfMatIsTransposed(A) ?
    bfMatColStride(A)/dtype_size :
    bfMatRowStride(A)/dtype_size;
}

static void
complexComplexMatMul(BfMat const *A, BfMat const *B, BfMat *C)
{
  assert(A->dtype != BF_DTYPE_MAT);

  assert(C->dtype == BF_DTYPE_COMPLEX);
  assert(!bfMatIsTransposed(C));

  BfComplex alpha = 1, beta = 0;

  enum CBLAS_TRANSPOSE A_trans = getCblasTranspose(A);
  enum CBLAS_TRANSPOSE B_trans = getCblasTranspose(B);

  BfSize m = bfMatNumRows(A);
  BfSize n = bfMatNumCols(B);
  BfSize k = bfMatNumCols(A);

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

  assert(bfMatIsAligned(A));
  assert(bfMatIsAligned(B));
  assert(bfMatIsAligned(C));

  if (bfMatIsTransposed(C))
    return BF_ERROR_INVALID_ARGUMENTS;

  if (bfMatNumRows(A) != bfMatNumRows(C))
    return BF_ERROR_INVALID_ARGUMENTS;

  if (bfMatNumCols(A) != bfMatNumRows(B))
    return BF_ERROR_INVALID_ARGUMENTS;

  if (bfMatNumCols(B) != bfMatNumCols(C))
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
  assert(A->dtype != BF_DTYPE_MAT);

  if (C->dtype != BF_DTYPE_COMPLEX)
    return BF_ERROR_INVALID_ARGUMENTS;

  BfReal scale;
  if (A->props & BF_MAT_PROP_DIAGONAL) {
    BfSize n = A->shape[1];
    BfVec row;
    for (BfSize i = 0; i < n; ++i) {
      bfCopyMatRow(B, i, C, i);
      row = bfGetMatRow(C, i);
      bfGetMatElt(A, i, i, &scale);
      scale = 1/scale;
      bfVecScale(&row, &scale);
    }
  }

  return BF_ERROR_NO_ERROR;
}

enum BfError
complexComplexMatSolve(BfMat const *A, BfMat const *B, BfMat *C)
{
  assert(A->dtype != BF_DTYPE_MAT);

  if (C->dtype != BF_DTYPE_COMPLEX)
    return BF_ERROR_INVALID_ARGUMENTS;

  if (A->props & BF_MAT_PROP_UNITARY) {
    BfSize m = A->shape[1], n = B->shape[1], k = A->shape[0];

    BfComplex alpha = 1, beta = 0;

    cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans, m, n, k,
                &alpha, A->data, m, B->data, n, &beta, C->data, n);

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

  bool A_trans = bfMatIsTransposed(A);
  bool B_trans = bfMatIsTransposed(B);
  bool C_trans = bfMatIsTransposed(C);

  /* check whether A and B have compatible shapes to form A\B */

  if ((!A_trans && !B_trans && A->shape[0] != B->shape[0]) ||
      ( A_trans && !B_trans && A->shape[1] != B->shape[0]) ||
      (!A_trans &&  B_trans && A->shape[0] != B->shape[1]) ||
      ( A_trans &&  B_trans && A->shape[1] != B->shape[1]))
    return BF_ERROR_INVALID_ARGUMENTS;

  /* check that A\B and C have the same shape */

  if ((!A_trans && !C_trans && A->shape[1] != C->shape[0]) ||
      ( A_trans && !C_trans && A->shape[0] != C->shape[0]) ||
      (!A_trans &&  C_trans && A->shape[1] != C->shape[1]) ||
      ( A_trans &&  C_trans && A->shape[0] != C->shape[1]))
    return BF_ERROR_INVALID_ARGUMENTS;

  if ((!B_trans && !C_trans && B->shape[1] != C->shape[1]) ||
      ( B_trans && !C_trans && B->shape[0] != C->shape[1]) ||
      (!B_trans &&  C_trans && B->shape[1] != C->shape[0]) ||
      ( B_trans &&  C_trans && B->shape[0] != C->shape[0]))
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

  enum BfError error = BF_ERROR_NO_ERROR;

  BfSize m = A->shape[0];
  BfSize n = A->shape[1];
  BfSize p = m < n ? m : n;

  /* init U */
  BfSize U_shape[] = {m, p};
  error = bfInitEmptyMat(U, BF_DTYPE_COMPLEX, BF_MAT_PROP_NONE, U_shape);
  if (error)
    return error;

  /* init S */
  BfSize S_shape[] = {p, p};
  error = bfInitEmptyMat(S, BF_DTYPE_REAL, BF_MAT_PROP_DIAGONAL, S_shape);
  if (error) {
    bfFreeMat(U);
    return error;
  }

  /* init Vt */
  BfSize Vt_shape[] = {p, n};
  error = bfInitEmptyMat(Vt, BF_DTYPE_COMPLEX, BF_MAT_PROP_NONE, Vt_shape);
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

  enum BfError error = BF_ERROR_NO_ERROR;

  lapack_int m = A->shape[0];
  lapack_int n = A->shape[1];

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

  enum BfError error;

  /* compute SVD of A */

  BfMat U, S, VH;

  error = bfInitEmptySvdMats(A, &U, &S, &VH);
  if (error)
    goto cleanup;

  error = bfComputeMatSvd(A, &U, &S, &VH);
  if (error)
    goto cleanup;

  /* compute tolerance and compute number of terms to
   * retain in pseudoinverse */

  BfSize m = A->shape[0], n = A->shape[1];
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

  BfMat WkH;

  error = bfInitEmptyMat(&WkH, UkH.dtype, BF_MAT_PROP_NONE, (BfSize[]) {k, n});
  if (error)
    goto cleanup;

  error = bfMatSolve(&Sk, &UkH, &WkH);
  if (error)
    goto cleanup;

  error = bfInitEmptyMat(pinv, A->dtype, A->props, (BfSize[]) {n, m});
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
  enum BfError error = BF_ERROR_NO_ERROR;

  BfReal const atol = BF_EPS_MACH;

  BfSize m = A->shape[0];
  BfSize n = A->shape[1];
  BfReal const rtol = (m > n ? m : n)*BF_EPS_MACH;

  BfMat A_pinv;
  error = bfComputePinv(A, atol, rtol, &A_pinv);
  if (error)
    goto cleanup;

  error = bfMatMul(&A_pinv, B, C);
  if (error)
    goto cleanup;

cleanup:
  bfFreeMat(&A_pinv);

  return error;
}
