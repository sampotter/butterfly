#include "mat.h"

#include <stdlib.h>
#include <string.h>

BfSize bfMatSize(BfMat const *mat)
{
  BfSize size = 1;
  for (BfSize i = 0; i < mat->ndim; ++i)
    size *= mat->shape[i];
  return size;
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

  mat->dtype = dtype;
  mat->ndim = ndim;

  memcpy(mat->shape, shape, sizeof(BfSize)*ndim);

  BfSize dtype_size;
  error |= bfSizeOfDtype(dtype, &dtype_size);
  if (error)
    return error;

  mat->data = malloc(dtype_size*bfMatSize(mat));
}
