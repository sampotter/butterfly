#include <bf/ptr_array.h>

#include <bf/assert.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mem.h>
#include <bf/util.h>

#include "macros.h"

/** Implementation: PtrArray */

BfPtrArray *bfPtrArrayGetView(BfPtrArray *ptrArray) {
  BF_ERROR_BEGIN();

  BfPtrArray *ptrArrayView = bfMemAlloc(1, sizeof(BfPtrArray));
  HANDLE_ERROR();

  *ptrArrayView = *ptrArray;

  ptrArrayView->isView = true;

  BF_ERROR_END() {
    BF_DIE();
  }

  return ptrArrayView;
}

BfPtrArray *bfPtrArrayCopy(BfPtrArray *ptrArray) {
  BF_ERROR_BEGIN();

  BfPtrArray *ptrArrayCopy = bfMemAlloc(1, sizeof(BfPtrArray));
  HANDLE_ERROR();

  ptrArrayCopy->data = bfMemAlloc(ptrArray->num_elts, sizeof(BfPtr));
  HANDLE_ERROR();

  bfMemCopy(ptrArray->data, ptrArray->num_elts, sizeof(BfPtr), ptrArrayCopy->data);

  ptrArrayCopy->capacity = ptrArray->num_elts;
  ptrArrayCopy->num_elts = ptrArray->num_elts;
  ptrArrayCopy->isView = false;

  BF_ERROR_END() {
    BF_DIE();
  }

  return ptrArrayCopy;
}

BfPtrArray *bfPtrArraySteal(BfPtrArray *ptrArray) {
  BF_ERROR_BEGIN();

  if (ptrArray->isView)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfPtrArray *ptrArrayStolen = bfMemAlloc(1, sizeof(BfPtrArray));
  HANDLE_ERROR();

  *ptrArrayStolen = *ptrArray;

  BF_ASSERT(!ptrArrayStolen->isView);
  ptrArray->isView = true;

  BF_ERROR_END() {
    BF_DIE();
  }

  return ptrArrayStolen;
}

BfPtrArray *bfPtrArrayNewWithDefaultCapacity() {
  BF_ERROR_BEGIN();

  BfPtrArray *ptrArray = bfMemAlloc(1, sizeof(BfPtrArray));
  HANDLE_ERROR();

  bfInitPtrArrayWithDefaultCapacity(ptrArray);
  HANDLE_ERROR();

  BF_ERROR_END() {
    ptrArray = NULL;
  }

  return ptrArray;
}

BfPtrArray bfGetUninitializedPtrArray() {
  return (BfPtrArray) {
    .data = NULL,
    .capacity = 0,
    .num_elts = 0,
    .isView = false
  };
}

void bfInitPtrArray(BfPtrArray *arr, BfSize capacity) {
  BF_ERROR_BEGIN();

  arr->data = bfMemAllocAndZero(capacity, sizeof(BfPtr));
  if (arr->data == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  arr->capacity = capacity;
  arr->num_elts = 0;
  arr->isView = false;

  BF_ERROR_END() {
    bfPtrArrayDeinit(arr);
  }
}

void bfInitPtrArrayWithDefaultCapacity(BfPtrArray *arr) {
  bfInitPtrArray(arr, BF_ARRAY_DEFAULT_CAPACITY);
}

void bfPtrArrayDeinit(BfPtrArray *arr) {
  bfMemFree(arr->data);
  bfMemZero(arr, 1, sizeof(BfPtrArray));
}

void bfPtrArrayDealloc(BfPtrArray **arr) {
  bfMemFree(*arr);
  *arr = NULL;
}

void bfPtrArrayDeinitAndDealloc(BfPtrArray **arr) {
  bfPtrArrayDeinit(*arr);
  bfPtrArrayDealloc(arr);
}

BfSize bfPtrArraySize(BfPtrArray const *arr) {
  return arr->num_elts;
}

bool bfPtrArrayIsEmpty(BfPtrArray const *arr) {
  return arr->num_elts == 0;
}

void extendPtrArray(BfPtrArray *arr, BfSize new_capacity) {
  BF_ERROR_BEGIN();

  BF_ASSERT(new_capacity > arr->capacity);

  void *new_data = bfMemRealloc(arr->data, new_capacity, sizeof(BfPtr));
  HANDLE_ERROR();

  arr->data = new_data;
  arr->capacity = new_capacity;

  BF_ERROR_END() {}
}

void bfPtrArrayAppend(BfPtrArray *arr, BfPtr ptr) {
  BF_ERROR_BEGIN();

  if (arr->num_elts == arr->capacity) {
    extendPtrArray(arr, 2*arr->capacity);
    HANDLE_ERROR();
  }

  arr->data[arr->num_elts++] = ptr;

  BF_ERROR_END() {}
}

BfPtr bfPtrArrayGet(BfPtrArray const *arr, BfSize pos) {
  BF_ERROR_BEGIN();

  BfPtr ptr = NULL;

  if (pos >= arr->num_elts)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  ptr = arr->data[pos];

  BF_ERROR_END() {}

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

BfPtrArray *bfPtrArrayGetRangeView(BfPtrArray const *arr, BfSize i0, BfSize i1) {
  BF_ERROR_BEGIN();

  if (i0 > i1)
    RAISE_ERROR(BF_ERROR_OUT_OF_RANGE);

  if (i1 > arr->num_elts)
    RAISE_ERROR(BF_ERROR_OUT_OF_RANGE);

  BfPtrArray *rangeView = bfPtrArrayGetView((BfPtrArray *)arr);
  HANDLE_ERROR();

  rangeView->data = &rangeView->data[i0];
  rangeView->num_elts = i1 - i0;

  /* The `capacity` variable should tell us how much space there is
   * until the end of the dynamically allocated array. We need to
   * shift it here to make sure it points to the right place. */
  rangeView->capacity -= i0;

  BF_ERROR_END() {
    BF_DIE();
  }

  return rangeView;
}

void bfMapPtrArray(BfPtrArray *arr, BfPtrFunc func, void *arg) {
  for (BfSize i = 0; i < arr->num_elts; ++i) {
    func(arr->data[i], arg);
    enum BfError error = bfGetError();
    if (error) {
      bfSetError(error);
      return;
    }
  }
}

BfPtr bfPtrArrayPopFirst(BfPtrArray *arr) {
  BF_ERROR_BEGIN();

  BfPtr ptr = NULL;

  if (arr->num_elts == 0)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  ptr = arr->data[0];
  bfMemMove(arr->data + 1, arr->num_elts - 1, sizeof(BfPtr), arr->data);
#if BF_DEBUG
  arr->data[arr->num_elts - 1] = NULL;
#endif
  --arr->num_elts;

  BF_ERROR_END() {
    BF_DIE();
  }

  return ptr;
}

BfPtr bfPtrArrayPopLast(BfPtrArray *arr) {
  BF_ERROR_BEGIN();

  BfPtr ptr = NULL;

  if (arr->num_elts == 0)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  ptr = arr->data[arr->num_elts - 1];
#if BF_DEBUG
  arr->data[arr->num_elts - 1] = NULL;
#endif
  --arr->num_elts;

  BF_ERROR_END() {
    BF_DIE();
  }

  return ptr;
}

typedef struct {
  BfPtrCmp ptrCmp;
} PtrCmpWrapper;

int ptrCmpWrapped(BfPtr ptr1, BfPtr ptr2, PtrCmpWrapper *wrapper) {
  return wrapper->ptrCmp(ptr1, ptr2);
}

void bfPtrArraySort(BfPtrArray *arr, BfPtrCmp ptrCmp, BfPtr aux) {
  (void)aux;
  PtrCmpWrapper wrapper = {.ptrCmp = ptrCmp};
  bfSort(arr->data, arr->num_elts, sizeof(BfPtr), (BfCompar)ptrCmpWrapped, &wrapper);
}

void bfPtrArrayReverse(BfPtrArray *arr) {
  BfSize n = bfPtrArraySize(arr);
  for (BfSize i = 0; i < n/2; ++i)
    SWAP(arr->data[i], arr->data[n - i - 1]);
}

void bfPtrArrayCopyData(BfPtrArray const *ptrArray, BfPtr *dst) {
  BF_ERROR_BEGIN();

  if (dst == NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  bfMemCopy(ptrArray->data, ptrArray->num_elts, sizeof(BfPtr), dst);

  BF_ERROR_END() {}
}

void bfPtrArrayRemove(BfPtrArray *ptrArray, BfSize i) {
  BF_ERROR_BEGIN();

  if (i >= bfPtrArraySize(ptrArray))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  bfMemMove(ptrArray->data + i + 1, ptrArray->num_elts - i - 1, sizeof(BfSize), ptrArray->data + i);

  --ptrArray->num_elts;

  BF_ERROR_END() {
    BF_DIE();
  }
}

/** Implementation: ConstPtrArray */

BfConstPtrArray *bfConstPtrArrayNewWithDefaultCapacity() {
  BF_ERROR_BEGIN();

  BfConstPtrArray *constPtrArray = bfMemAlloc(1, sizeof(BfConstPtrArray));
  if (constPtrArray == NULL)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  bfConstPtrArrayInitWithDefaultCapacity(constPtrArray);
  HANDLE_ERROR();

  BF_ERROR_END() {
    constPtrArray = NULL;
  }

  return constPtrArray;
}

void bfConstPtrArrayInit(BfConstPtrArray *arr, BfSize capacity) {
  BF_ERROR_BEGIN();

  arr->data = bfMemAllocAndZero(capacity, sizeof(BfPtr));
  if (arr->data == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  arr->capacity = capacity;
  arr->num_elts = 0;
  arr->isView = true;

  BF_ERROR_END() {
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
  BF_ASSERT(new_capacity > arr->capacity);

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
  BF_ERROR_BEGIN();

  BfConstPtr ptr = NULL;

  if (pos >= arr->num_elts)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  ptr = arr->data[pos];

  BF_ERROR_END() {}

  return ptr;
}

BfConstPtr bfConstPtrArrayPopLast(BfConstPtrArray *arr) {
  BfConstPtr ptr = arr->data[--arr->num_elts];
  arr->data[arr->num_elts] = NULL;
  return ptr;
}

BfConstPtr bfConstPtrArrayGetFirst(BfConstPtrArray const *arr) {
  BF_ERROR_BEGIN();

  BfConstPtr constPtr = NULL;

  if (bfConstPtrArrayIsEmpty(arr))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  constPtr = arr->data[0];

  BF_ERROR_END() {}

  return constPtr;
}

BfConstPtr bfConstPtrArrayGetLast(BfConstPtrArray const *arr) {
  BF_ERROR_BEGIN();

  BfConstPtr constPtr = NULL;

  if (bfConstPtrArrayIsEmpty(arr))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  constPtr = arr->data[arr->num_elts - 1];

  BF_ERROR_END() {}

  return constPtr;
}

void bfConstPtrArrayExtend(BfConstPtrArray *arr, BfConstPtrArray const *otherArr) {
  BF_ERROR_BEGIN();

  for (BfSize i = 0; i < bfConstPtrArraySize(otherArr); ++i) {
    bfConstPtrArrayAppend(arr, bfConstPtrArrayGet(otherArr, i));
    HANDLE_ERROR();
  }

  BF_ERROR_END() {
    BF_DIE(); // roll back changes?
  }
}
