#include "mat.h"

#include <stdlib.h>

BfSize bfMatSize(BfMat const *mat)
{
  BfSize size = 1;
  for (BfSize i = 0; i < mat->ndim; ++i)
    size *= mat->shape[i];
  return size;
}

enum BfError bfMatNumBytes(BfMat const *mat, BfSize *num_bytes)
{
  enum BfError error = BF_ERROR_NO_ERROR;

  BfSize dtype_size;
  error |= bfSizeOfDtype(mat->dtype, &dtype_size);
  if (error)
    return error;

  *num_bytes = dtype_size*bfMatSize(mat);

  return error;
}

enum BfError bfFreeMat(BfMat *mat)
{
  free(mat->data);
}

enum BfError
bfMakeEmptyMat(BfMat *mat,
               enum BfDtypes dtype, BfSize ndim, BfSize const *shape)
{
  enum BfError error = BF_ERROR_NO_ERROR;

  BfSize dtype_size;
  error |= bfSizeOfDtype(dtype, &dtype_size);
  if (error)
    return error;

  mat->dtype = dtype;
  mat->ndim = ndim;

  mat->shape = malloc(sizeof(BfSize)*ndim);
  for (BfSize i = 0; i < ndim; ++i)
    mat->shape[i] = shape[i];

  mat->data = malloc(dtype_size*bfMatSize(mat));
}
