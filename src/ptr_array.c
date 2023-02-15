#include <bf/ptr_array.h>

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <bf/error.h>
#include <bf/error_macros.h>

#include "macros.h"

/** Implementation: PtrArray */

BfPtrArray bfGetUninitializedPtrArray() {
  return (BfPtrArray) {
    .flags = BF_PTR_ARRAY_FLAG_NONE,
    .data = NULL,
    .capacity = 0,
    .num_elts = 0
  };
}

void bfInitPtrArray(BfPtrArray *arr, BfSize capacity) {
  BEGIN_ERROR_HANDLING();

  arr->flags = BF_PTR_ARRAY_FLAG_NONE;

  arr->data = malloc(capacity*sizeof(BfPtr));
  if (arr->data == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  arr->capacity = capacity;
  arr->num_elts = 0;

  END_ERROR_HANDLING() {
    bfPtrArrayDeinit(arr);
  }
}

void bfInitPtrArrayWithDefaultCapacity(BfPtrArray *arr) {
  bfInitPtrArray(arr, BF_ARRAY_DEFAULT_CAPACITY);
}

void bfMakeEmptyPtrArrayView(BfPtrArray *arr) {
  arr->flags = BF_PTR_ARRAY_FLAG_NONE;
  arr->data = NULL;
  arr->capacity = 0;
  arr->num_elts = 0;
}

void bfPtrArrayDeinit(BfPtrArray *arr) {
  free(arr->data);
  memset(arr, 0x0, sizeof(BfPtrArray));
}

BfPtrArray bfPtrArrayCopy(BfPtrArray *arr) {
  BEGIN_ERROR_HANDLING();

  BfPtrArray arrCopy = {
    .flags = BF_PTR_ARRAY_FLAG_NONE,
    .data = malloc(arr->capacity*sizeof(BfPtr)),
    .capacity = arr->capacity,
    .num_elts = arr->num_elts
  };

  if (arrCopy.data == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  memcpy(arrCopy.data, arr->data, arr->num_elts*sizeof(BfPtr));

  END_ERROR_HANDLING()
    free(arrCopy.data);

  return arrCopy;
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

void bfPtrArrayAppend(BfPtrArray *arr, BfPtr ptr) {
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

BfPtr bfPtrArrayGet(BfPtrArray const *arr, BfSize pos) {
  BEGIN_ERROR_HANDLING();

  BfPtr ptr = NULL;

  if (pos >= arr->num_elts)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  ptr = arr->data[pos];

  END_ERROR_HANDLING() {}

  return ptr;
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

BfPtr bfPtrArrayPopLast(BfPtrArray *arr) {
  BfPtr ptr = arr->data[--arr->num_elts];
  arr->data[arr->num_elts] = NULL;
  return ptr;
}

void bfPtrArraySort(BfPtrArray *arr, BfPtrCmp ptrCmp) {
  qsort(arr->data, arr->num_elts, sizeof(BfPtr), (__compar_fn_t)ptrCmp);
}

void bfPtrArrayReverse(BfPtrArray *arr) {
  BfSize n = bfPtrArraySize(arr);
  for (BfSize i = 0; i < n/2; ++i)
    SWAP(arr->data[i], arr->data[n - i - 1]);
}

void bfPtrArrayCopyData(BfPtrArray const *ptrArray, BfPtr *dst) {
  BEGIN_ERROR_HANDLING();

  if (dst == NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  memcpy(dst, ptrArray->data, ptrArray->num_elts*sizeof(BfPtr));

  END_ERROR_HANDLING() {}
}

/** Implementation: ConstPtrArray */

void bfConstPtrArrayDeinit(BfConstPtrArray *arr) {
  free(arr->data);
  memset(arr, 0x0, sizeof(BfConstPtrArray));
}

bool bfConstPtrArrayIsEmpty(BfConstPtrArray const *arr) {
  return arr->num_elts == 0;
}

void extendConstPtrArray(BfConstPtrArray *arr, BfSize new_capacity) {
  assert(new_capacity > arr->capacity);

  void *new_data = realloc(arr->data, sizeof(BfConstPtr)*new_capacity);
  if (new_data == NULL) {
    bfSetError(BF_ERROR_MEMORY_ERROR);
    return;
  }

  arr->data = new_data;
  arr->capacity = new_capacity;
}

void bfConstPtrArrayAppend(BfConstPtrArray *arr, BfConstPtr constPtr) {
  if (arr->num_elts == arr->capacity) {
    extendConstPtrArray(arr, 2*arr->capacity);
    enum BfError error = bfGetError();
    if (error) {
      bfSetError(error);
      return;
    }
  }
  arr->data[arr->num_elts++] = constPtr;
}

BfConstPtr bfConstPtrArrayPopLast(BfConstPtrArray *arr) {
  BfConstPtr ptr = arr->data[--arr->num_elts];
  arr->data[arr->num_elts] = NULL;
  return ptr;
}
