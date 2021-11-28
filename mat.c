#include "mat.h"

#include <stdlib.h>

#include <lapack.h>

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

  lapack_int m = A->shape[0], n = A->shape[1];
  lapack_int lwork = -1, info;

  BfComplex wkopt, *work;
  BfReal *rwork = malloc(5*(m < n ? m : n)*sizeof(BfReal));

  /* calculate the optimal workspace */
  LAPACK_zgesvd(
    "S", "S", &m, &n, (BfComplex *)A->data, &m, (BfReal *)S->data,
    (BfComplex *)U->data, &m, (BfComplex *)Vt->data, &n, &wkopt,
    &lwork, rwork, &info);

  lwork = (int)creal(wkopt);
  work = malloc(lwork*sizeof(BfComplex));

  /* compute the SVD */
  LAPACK_zgesvd(
    "S", "S", &m, &n, (BfComplex *)A->data, &m, (BfReal *)S->data,
    (BfComplex *)U->data, &m, (BfComplex *)Vt->data, &n, work,
    &lwork, rwork, &info);

  /* check for errors */
  if (info > 0)
    error |= BF_ERROR_RUNTIME_ERROR;

  free(work);

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
