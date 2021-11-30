#include "mat.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <cblas.h>
#include <lapacke.h>

#include "rand.h"

BfSize bfMatSize(BfMat const *A)
{
  return A->shape[0]*A->shape[1];
}

enum BfError bfMatNumBytes(BfMat const *A, BfSize *nbytes)
{
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

enum BfError bfFreeMat(BfMat *A)
{
  free(A->data);

  return BF_ERROR_NO_ERROR;
}

enum BfError
bfInitEmptyMat(BfMat *A, enum BfDtypes dtype, enum BfMatProps props,
               BfSize const shape[2])
{
  enum BfError error = BF_ERROR_NO_ERROR;

  A->dtype = dtype;
  A->props = props;
  A->shape[0] = shape[0];
  A->shape[1] = shape[1];

  BfSize nbytes;
  error = bfMatNumBytes(A, &nbytes);
  if (error)
    return error;

  A->data = malloc(nbytes);

  return error;
}

enum BfError
bfSaveMat(BfMat const *A, char const *path)
{
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

enum BfError
bfGetMatRow(BfMat const *A, BfSize i, void **data)
{
  BfSize m = A->shape[0], n = A->shape[1];

  if (i >= m)
    return BF_ERROR_INVALID_ARGUMENTS;

  BfSize nbytes;
  enum BfError error = bfSizeOfDtype(A->dtype, &nbytes);
  if (error)
    return error;

  BfSize row_nbytes = nbytes*n;

  *data = A->data + row_nbytes*i;

  return BF_ERROR_NO_ERROR;
}

enum BfError
bfSetMatRow(BfMat *A, BfSize i, void const *data)
{
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
complexComplexMatMul(BfMat const *A, BfMat const *B, BfMat *C)
{
  if (C->dtype != BF_DTYPE_COMPLEX)
    return BF_ERROR_INVALID_ARGUMENTS;

  if (A->props & BF_MAT_PROP_CONJ_TRANS)
    return BF_ERROR_NOT_IMPLEMENTED;

  if (B->props & BF_MAT_PROP_CONJ_TRANS)
    return BF_ERROR_NOT_IMPLEMENTED;

  BfSize m = A->shape[0], n = B->shape[1], k = A->shape[1];

  BfComplex alpha = 1, beta = 0;

  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
              &alpha, A->data, k, B->data, n, &beta, C->data, n);

  return BF_ERROR_NO_ERROR;
}

enum BfError
bfMatMul(BfMat const *A, BfMat const *B, BfMat *C)
{
  if (A->shape[0] != C->shape[0])
    return BF_ERROR_INVALID_ARGUMENTS;

  if (A->shape[1] != B->shape[0])
    return BF_ERROR_INVALID_ARGUMENTS;

  if (B->shape[1] != C->shape[1])
    return BF_ERROR_INVALID_ARGUMENTS;

  if (A->dtype == BF_DTYPE_COMPLEX && B->dtype == BF_DTYPE_COMPLEX)
    return complexComplexMatMul(A, B, C);
  else
    return BF_ERROR_NOT_IMPLEMENTED;

  return BF_ERROR_NO_ERROR;
}

enum BfError
realComplexMatSolve(BfMat const *A, BfMat const *B, BfMat *C)
{
  if (C->dtype != BF_DTYPE_COMPLEX)
    return BF_ERROR_INVALID_ARGUMENTS;

  if (A->props & BF_MAT_PROP_DIAGONAL) {
    BfSize m = A->shape[0], n = B->shape[1], k = A->shape[1];

    if (m != k)
      return BF_ERROR_NOT_IMPLEMENTED;

    BfReal a;
    BfComplex *row;
    for (BfSize i = 0; i < m; ++i) {
      bfGetMatRow(B, i, (void **)&row);
      bfSetMatRow(C, i, row);
      bfGetMatRow(C, i, (void **)&row);
      a = *((BfReal const *)A->data + i);
      cblas_zdscal(n, 1/a, row, 1);
    }

    return BF_ERROR_NO_ERROR;
  }

  return BF_ERROR_NOT_IMPLEMENTED;
}

enum BfError
complexComplexMatSolve(BfMat const *A, BfMat const *B, BfMat *C)
{
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
  if (A->shape[0] != B->shape[0])
    return BF_ERROR_INVALID_ARGUMENTS;

  if (A->shape[1] != C->shape[0])
    return BF_ERROR_INVALID_ARGUMENTS;

  if (B->shape[1] != C->shape[1])
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
  switch (A->dtype) {
  case BF_DTYPE_COMPLEX:
    return computeMatSvdComplex(A, U, S, Vt);
  default:
    return BF_ERROR_NOT_IMPLEMENTED;
  }
}
