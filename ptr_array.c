#include "ptr_array.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "error_macros.h"

BfPtrArray bfGetUninitializedPtrArray() {
  return (BfPtrArray) {
    .flags = BF_PTR_ARRAY_FLAG_NONE,
    .data = NULL,
    .capacity = 0,
    .num_elts = 0
  };
}

void
bfInitPtrArray(BfPtrArray *arr, BfSize capacity)
{
  arr->flags = BF_PTR_ARRAY_FLAG_NONE;

  arr->data = malloc(capacity*sizeof(BfPtr));
  if (arr->data == NULL) {
    bfSetError(BF_ERROR_MEMORY_ERROR);
    return;
  }

  arr->capacity = capacity;
  arr->num_elts = 0;
}

void
bfInitPtrArrayWithDefaultCapacity(BfPtrArray *arr)
{
  bfInitPtrArray(arr, BF_ARRAY_DEFAULT_CAPACITY);
}

void bfMakeEmptyPtrArrayView(BfPtrArray *arr)
{
  arr->flags = BF_PTR_ARRAY_FLAG_NONE;
  arr->data = NULL;
  arr->capacity = 0;
  arr->num_elts = 0;
}

void
bfFreePtrArray(BfPtrArray *arr)
{
  free(arr->data);
}

BfSize bfPtrArraySize(BfPtrArray const *arr) {
  return arr->num_elts;
}

bool bfPtrArrayIsEmpty(BfPtrArray const *arr) {
  return arr->num_elts == 0;
}

void extendPtrArray(BfPtrArray *arr, BfSize new_capacity) {
  assert(new_capacity > arr->capacity);

  void *new_data = realloc(arr->data, sizeof(BfPtr)*new_capacity);
  if (new_data == NULL) {
    bfSetError(BF_ERROR_MEMORY_ERROR);
    return;
  }

  arr->data = new_data;
  arr->capacity = new_capacity;
}

void
bfPtrArrayAppend(BfPtrArray *arr, BfPtr ptr)
{
  if (arr->num_elts == arr->capacity) {
    extendPtrArray(arr, 2*arr->capacity);
    enum BfError error = bfGetError();
    if (error) {
      bfSetError(error);
      return;
    }
  }

  arr->data[arr->num_elts++] = ptr;
}

void
bfPtrArrayGet(BfPtrArray const *arr, BfSize pos, BfPtr *elt)
{
  if (pos >= arr->num_elts) {
    bfSetError(BF_ERROR_INVALID_ARGUMENTS);
    return;
  }

  *elt = arr->data[pos];
}

void
bfPtrArrayGetFirst(BfPtrArray const *arr, BfPtr *ptr)
{
  if (arr->num_elts == 0) {
    bfSetError(BF_ERROR_INVALID_ARGUMENTS);
    return;
  }

  *ptr = arr->data[0];
}

void
bfPtrArrayGetLast(BfPtrArray const *arr, BfPtr *ptr)
{
  if (arr->num_elts == 0) {
    bfSetError(BF_ERROR_INVALID_ARGUMENTS);
    return;
  }

  *ptr = arr->data[arr->num_elts - 1];
}

void
bfPtrArrayGetRangeView(BfPtrArray const *arr, BfSize start, BfSize end,
                       BfPtrArray *view)
{
  if (start > end) {
    bfSetError(BF_ERROR_INVALID_ARGUMENTS);
    return;
  }

  if (end > arr->num_elts) {
    bfSetError(BF_ERROR_INVALID_ARGUMENTS);
    return;
  }

  view->flags = BF_PTR_ARRAY_FLAG_VIEW;
  view->data = arr->data + start;
  view->capacity = end - start;
  view->num_elts = end - start;
}

void
bfMapPtrArray(BfPtrArray *arr, BfPtrFunc func, void *arg)
{
  for (BfSize i = 0; i < arr->num_elts; ++i) {
    func(arr->data[i], arg);
    enum BfError error = bfGetError();
    if (error) {
      bfSetError(error);
      return;
    }
  }
}
