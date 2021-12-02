#include "ptr_array.h"

#include <stdlib.h>
#include <string.h>

enum BfError
bfInitPtrArray(BfPtrArray *arr, BfSize capacity)
{
  arr->flags = BF_PTR_ARRAY_FLAG_NONE;

  arr->data = malloc(capacity*sizeof(BfPtr));
  if (arr->data == NULL)
    return BF_ERROR_MEMORY_ERROR;

  arr->capacity = capacity;
  arr->num_elts = 0;

  return BF_ERROR_NO_ERROR;
}

enum BfError
bfInitPtrArrayWithDefaultCapacity(BfPtrArray *arr)
{
  return bfInitPtrArray(arr, BF_ARRAY_DEFAULT_CAPACITY);
}

void bfMakeEmptyPtrArrayView(BfPtrArray *arr)
{
  arr->flags = BF_PTR_ARRAY_FLAG_NONE;
  arr->data = NULL;
  arr->capacity = 0;
  arr->num_elts = 0;
}

enum BfError
bfFreePtrArray(BfPtrArray *arr)
{
  free(arr->data);

  return BF_ERROR_NO_ERROR;
}

BfSize bfPtrArraySize(BfPtrArray const *arr) {
  return arr->num_elts;
}

bool bfPtrArrayIsEmpty(BfPtrArray const *arr) {
  return arr->num_elts == 0;
}

enum BfError extendPtrArray(BfPtrArray *arr, BfSize new_capacity) {
  if (new_capacity <= arr->capacity)
    return BF_ERROR_INVALID_ARGUMENTS;

  void *new_data = reallocarray(arr->data, new_capacity, sizeof(BfPtr));
  if (new_data == NULL)
    return BF_ERROR_MEMORY_ERROR;

  arr->data = new_data;
  arr->capacity = new_capacity;

  return BF_ERROR_NO_ERROR;
}

enum BfError
bfPtrArrayAppend(BfPtrArray *arr, BfPtr ptr)
{
  enum BfError error = BF_ERROR_NO_ERROR;

  if (arr->num_elts == arr->capacity) {
    error = extendPtrArray(arr, 2*arr->capacity);
    if (error)
      return error;
  }

  arr->data[arr->num_elts++] = ptr;

  return error;
}

enum BfError
bfPtrArrayGet(BfPtrArray const *arr, BfSize pos, BfPtr *elt)
{
  if (pos >= arr->num_elts)
    return BF_ERROR_INVALID_ARGUMENTS;

  *elt = arr->data[pos];

  return BF_ERROR_NO_ERROR;
}

enum BfError
bfPtrArrayGetFirst(BfPtrArray const *arr, BfPtr *ptr)
{
  if (arr->num_elts == 0)
    return BF_ERROR_INVALID_ARGUMENTS;

  *ptr = arr->data[0];

  return BF_ERROR_NO_ERROR;
}

enum BfError
bfPtrArrayGetLast(BfPtrArray const *arr, BfPtr *ptr)
{
  if (arr->num_elts == 0)
    return BF_ERROR_INVALID_ARGUMENTS;

  *ptr = arr->data[arr->num_elts - 1];

  return BF_ERROR_NO_ERROR;
}

enum BfError
bfPtrArrayGetRangeView(BfPtrArray const *arr, BfSize start, BfSize end,
                       BfPtrArray *view)
{
  if (start > end)
    return BF_ERROR_INVALID_ARGUMENTS;

  if (end > arr->num_elts)
    return BF_ERROR_INVALID_ARGUMENTS;

  view->flags = BF_PTR_ARRAY_FLAG_VIEW;
  view->data = arr->data + start;
  view->capacity = end - start;
  view->num_elts = end - start;

  return BF_ERROR_NO_ERROR;
}

enum BfError
bfMapPtrArray(BfPtrArray *arr, BfPtrFunc func, void *arg)
{
  enum BfError error = BF_ERROR_NO_ERROR;

  for (BfSize i = 0; i < arr->num_elts; ++i) {
    error = func(arr->data[i], arg);
    if (error)
      break;
  }

  return error;
}
