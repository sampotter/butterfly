#include <bf/ptr_array.h>

#include <assert.h>

#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mem.h>
#include <bf/util.h>

#include "macros.h"

/** Implementation: PtrArray */

BfPtrArray *bfPtrArrayNewWithDefaultCapacity() {
  BEGIN_ERROR_HANDLING();

  BfPtrArray *ptrArray = bfMemAlloc(1, sizeof(BfPtrArray));
  if (ptrArray == NULL)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  bfInitPtrArrayWithDefaultCapacity(ptrArray);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {
    ptrArray = NULL;
  }

  return ptrArray;
}

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

  arr->data = bfMemAllocAndZero(capacity, sizeof(BfPtr));
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
  bfMemFree(arr->data);
  bfMemZero(arr, 1, sizeof(BfPtrArray));
}

void bfPtrArrayDealloc(BfPtrArray **arr) {
  bfMemFree(*arr);
  *arr = NULL;
}

void bfPtrArrayDelete(BfPtrArray **arr) {
  bfPtrArrayDeinit(*arr);
  bfPtrArrayDealloc(arr);
}

BfPtrArray bfPtrArrayCopy(BfPtrArray *arr) {
  BEGIN_ERROR_HANDLING();

  BfPtrArray arrCopy = {
    .flags = BF_PTR_ARRAY_FLAG_NONE,
    .data = bfMemAllocAndZero(arr->capacity, sizeof(BfPtr)),
    .capacity = arr->capacity,
    .num_elts = arr->num_elts
  };

  if (arrCopy.data == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  bfMemCopy(arr->data, arr->num_elts, sizeof(BfPtr), arrCopy.data);

  END_ERROR_HANDLING()
    bfMemFree(arrCopy.data);

  return arrCopy;
}

BfSize bfPtrArraySize(BfPtrArray const *arr) {
  return arr->num_elts;
}

bool bfPtrArrayIsEmpty(BfPtrArray const *arr) {
  return arr->num_elts == 0;
}

void extendPtrArray(BfPtrArray *arr, BfSize new_capacity) {
  BEGIN_ERROR_HANDLING();

  assert(new_capacity > arr->capacity);

  void *new_data = bfMemRealloc(arr->data, new_capacity, sizeof(BfPtr));
  HANDLE_ERROR();

  arr->data = new_data;
  arr->capacity = new_capacity;

  END_ERROR_HANDLING() {}
}

void bfPtrArrayAppend(BfPtrArray *arr, BfPtr ptr) {
  BEGIN_ERROR_HANDLING();

  if (arr->num_elts == arr->capacity) {
    extendPtrArray(arr, 2*arr->capacity);
    HANDLE_ERROR();
  }

  arr->data[arr->num_elts++] = ptr;

  END_ERROR_HANDLING() {}
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

void bfPtrArraySort(BfPtrArray *arr, BfPtrCmp ptrCmp, BfPtr aux) {
  bfSort(arr->data, arr->num_elts, sizeof(BfPtr), (BfCompar)ptrCmp, aux);
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

  bfMemCopy(ptrArray->data, ptrArray->num_elts, sizeof(BfPtr), dst);

  END_ERROR_HANDLING() {}
}

/** Implementation: ConstPtrArray */

BfConstPtrArray *bfConstPtrArrayNewWithDefaultCapacity() {
  BEGIN_ERROR_HANDLING();

  BfConstPtrArray *constPtrArray = bfMemAlloc(1, sizeof(BfConstPtrArray));
  if (constPtrArray == NULL)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  bfConstPtrArrayInitWithDefaultCapacity(constPtrArray);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {
    constPtrArray = NULL;
  }

  return constPtrArray;
}

void bfConstPtrArrayInit(BfConstPtrArray *arr, BfSize capacity) {
  BEGIN_ERROR_HANDLING();

  arr->flags = BF_PTR_ARRAY_FLAG_NONE;

  arr->data = bfMemAllocAndZero(capacity, sizeof(BfPtr));
  if (arr->data == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  arr->capacity = capacity;
  arr->num_elts = 0;

  END_ERROR_HANDLING() {
    bfConstPtrArrayDeinit(arr);
  }
}

void bfConstPtrArrayInitWithDefaultCapacity(BfConstPtrArray *arr) {
  bfConstPtrArrayInit(arr, BF_ARRAY_DEFAULT_CAPACITY);
}

void bfConstPtrArrayDeinit(BfConstPtrArray *arr) {
  bfMemFree(arr->data);
  bfMemZero(arr, 1, sizeof(BfConstPtrArray));
}

void bfConstPtrArrayDealloc(BfConstPtrArray **arr) {
  bfMemFree(*arr);
  *arr = NULL;
}

void bfConstPtrArrayDeinitAndDealloc(BfConstPtrArray **arr) {
  bfConstPtrArrayDeinit(*arr);
  bfConstPtrArrayDealloc(arr);
}

bool bfConstPtrArrayIsEmpty(BfConstPtrArray const *arr) {
  return arr->num_elts == 0;
}

BfSize bfConstPtrArraySize(BfConstPtrArray const *arr) {
  return arr->num_elts;
}

void extendConstPtrArray(BfConstPtrArray *arr, BfSize new_capacity) {
  assert(new_capacity > arr->capacity);

  void *new_data = bfMemRealloc(arr->data, new_capacity, sizeof(BfConstPtr));
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

BfConstPtr bfConstPtrArrayGet(BfConstPtrArray const *arr, BfSize pos) {
  BEGIN_ERROR_HANDLING();

  BfConstPtr ptr = NULL;

  if (pos >= arr->num_elts)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  ptr = arr->data[pos];

  END_ERROR_HANDLING() {}

  return ptr;
}

BfConstPtr bfConstPtrArrayPopLast(BfConstPtrArray *arr) {
  BfConstPtr ptr = arr->data[--arr->num_elts];
  arr->data[arr->num_elts] = NULL;
  return ptr;
}

BfConstPtr bfConstPtrArrayGetFirst(BfConstPtrArray const *arr) {
  BEGIN_ERROR_HANDLING();

  BfConstPtr constPtr = NULL;

  if (bfConstPtrArrayIsEmpty(arr))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  constPtr = arr->data[0];

  END_ERROR_HANDLING() {}

  return constPtr;
}

BfConstPtr bfConstPtrArrayGetLast(BfConstPtrArray const *arr) {
  BEGIN_ERROR_HANDLING();

  BfConstPtr constPtr = NULL;

  if (bfConstPtrArrayIsEmpty(arr))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  constPtr = arr->data[arr->num_elts - 1];

  END_ERROR_HANDLING() {}

  return constPtr;
}

void bfConstPtrArrayExtend(BfConstPtrArray *arr, BfConstPtrArray const *otherArr) {
  BEGIN_ERROR_HANDLING();

  for (BfSize i = 0; i < bfConstPtrArraySize(otherArr); ++i) {
    bfConstPtrArrayAppend(arr, bfConstPtrArrayGet(otherArr, i));
    HANDLE_ERROR();
  }

  END_ERROR_HANDLING() {
    assert(false); // roll back changes?
  }
}
