#include <bf/array.h>

#include <bf/assert.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mem.h>

#include <stdio.h>

#include "macros.h"

#define ELT_PTR(arr, i) arr->data + (i)*arr->eltSize

struct BfArray {
  BfSize eltSize;
  BfSize capacity;
  BfSize size;
  BfByte *data;
  bool isView;
};

static void invalidate(BfArray *array) {
  array->eltSize = BF_SIZE_BAD_VALUE;
  array->capacity = BF_SIZE_BAD_VALUE;
  array->size = BF_SIZE_BAD_VALUE;
  array->data = NULL;
  array->isView = false;
}

BfArray *bfArrayCopy(BfArray const *array) {
  BF_ERROR_BEGIN();

  BfArray *arrayCopy = bfArrayNewUninitialized();
  HANDLE_ERROR();

  arrayCopy->eltSize = array->eltSize;
  arrayCopy->capacity = array->capacity;
  arrayCopy->size = array->size;

  arrayCopy->data = bfMemAllocCopy(array->data, array->size, array->eltSize);
  HANDLE_ERROR();

  arrayCopy->isView = false;

  BF_ERROR_END() {
    BF_DIE();
  }

  return arrayCopy;
}

BfArray *bfArrayNewUninitialized(void) {
  BF_ERROR_BEGIN();

  BfArray *array = bfMemAlloc(1, sizeof(BfArray));
  HANDLE_ERROR();

  invalidate(array);

  BF_ERROR_END() {
    BF_DIE();
  }

  return array;
}

BfArray *bfArrayNewEmpty(BfSize eltSize) {
  BF_ERROR_BEGIN();

  BfArray *array = bfArrayNewUninitialized();
  HANDLE_ERROR();

  bfArrayInitEmpty(array, eltSize);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  return array;
}

BfArray *bfArrayNewWithValue(BfSize eltSize, BfSize size, BfConstPtr eltPtr) {
  BF_ERROR_BEGIN();

  BfArray *array = bfArrayNewUninitialized();
  HANDLE_ERROR();

  bfArrayInitWithValue(array, eltSize, size, eltPtr);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  return array;
}

void bfArrayInitEmpty(BfArray *array, BfSize eltSize) {
  BF_ERROR_BEGIN();

  array->eltSize = eltSize;
  array->capacity = BF_ARRAY_DEFAULT_CAPACITY;
  array->size = 0;

  array->data = bfMemAlloc(array->capacity, eltSize);
  HANDLE_ERROR();

  array->isView = false;

  BF_ERROR_END() {
    BF_DIE();
  }
}

void bfArrayInitWithValue(BfArray *array, BfSize eltSize, BfSize size, BfConstPtr eltPtr) {
  BF_ERROR_BEGIN();

  array->eltSize = eltSize;
  array->capacity = size;
  array->size = size;

  array->data = bfMemAlloc(array->capacity, array->eltSize);
  HANDLE_ERROR();

  array->isView = false;

  for (BfSize i = 0; i < array->size; ++i)
    bfMemCopy(eltPtr, 1, array->eltSize, ELT_PTR(array, i));

  BF_ERROR_END() {
    BF_DIE();
  }
}

void bfArrayDeinit(BfArray *array) {
  bfMemFree(array->data);

  invalidate(array);
}

void bfArrayDealloc(BfArray **array) {
  bfMemFree(*array);
  *array = NULL;
}

void bfArrayDeinitAndDealloc(BfArray **array) {
  bfArrayDeinit(*array);
  bfArrayDealloc(array);
}

BfSize bfArrayGetSize(BfArray const *array) {
  return array->size;
}

void bfArrayGet(BfArray const *array, BfSize i, BfPtr eltPtr) {
  BF_ERROR_BEGIN();

  if (i >= array->size)
    RAISE_ERROR(BF_ERROR_OUT_OF_RANGE);

  bfMemCopy(ELT_PTR(array, i), 1, array->eltSize, eltPtr);

  BF_ERROR_END() {
    BF_DIE();
  }
}

BfPtr bfArrayGetPtr(BfArray *array, BfSize i) {
  BF_ERROR_BEGIN();

  BfPtr ptr = NULL;

  if (i >= array->size)
    RAISE_ERROR(BF_ERROR_OUT_OF_RANGE);

  ptr = ELT_PTR(array, i);

  BF_ERROR_END() {
    BF_DIE();
  }

  return ptr;
}

BfSize bfArrayFindSorted(BfArray const *array, BfConstPtr eltPtr, BfCompar compar) {
  BF_ERROR_BEGIN();

  BfSize n = array->size;

  if (n == 0)
    return 0;

  if (n == 1)
    return compar(eltPtr, ELT_PTR(array, 0), NULL) <= 0 ? 0 : 1;

  BfSize i0 = 0;
  BfSize i1 = n;

#if BF_DEBUG
  BfSize iterCount = 0;
#endif

  BfSize i = BF_SIZE_BAD_VALUE;

  while (true) {
#if BF_DEBUG
    BF_ASSERT(iterCount <= n);
#endif
    BF_ASSERT(i0 <= i1 && i1 <= n);

    i = (i0 + i1)/2;
    BF_ASSERT(i0 <= i && i <= i1);

    /* If the index is zero, we can safely exit: */
    if (i == 0 || i == n) break;
    BF_ASSERT(0 < i && i < n);

    /* Evaluate the comparator between the `i-1`st element and the
     * target element. If `cmpPrev < 0`, then `arr[i-1] < tgt`; if
     * `cmpPrev > 0`, then `arr[i-1] > tgt`, and if `cmpPrev == 0`,
     * then `arr[i-1] == tgt`. */
    int cmpPrev = compar(ELT_PTR(array, i - 1), eltPtr, NULL);
    HANDLE_ERROR();

    /* Similar to `cmpPrev`, but for the `i`th element instead of the
     * `i-1`st. */
    int cmp = compar(ELT_PTR(array, i), eltPtr, NULL);
    HANDLE_ERROR();

    /* If `arr[i-1] < tgt` (i.e., `cmpPrev < 0`) and if `tgt <=
     * arr[i]` (i.e., `cmp >= 0`), then `i` is the lower bound: */
    if (cmpPrev < 0 && cmp >= 0) break;

    /* If `tgt <= arr[i - 1]` (`cmpPrev >= 0`), then set the bracket
     * to `[i0, i)`: */
    if (cmpPrev >= 0)
      i1 = i;

    /* If `tgt > arr[i]` (`cmp < 0`), then we update the bracket to
     * `[i0 + 1, i1)`: */
    if (cmp < 0)
      i0 = i + 1;

#if BF_DEBUG
    ++iterCount;
#endif
  }

  BF_ERROR_END() {
    BF_DIE();
  }

  BF_ASSERT(!BF_SIZE_OK(i) || i <= n);

  return i;
}

static void expand(BfArray *array) {
  BF_ERROR_BEGIN();

  BfSize newCapacity = 2*array->capacity;

  BfPtr newData = bfMemRealloc(array->data, newCapacity, array->eltSize);
  HANDLE_ERROR();

  array->capacity = newCapacity;
  array->data = newData;

  BF_ERROR_END() {
    BF_DIE();
  }
}

void bfArrayInsert(BfArray *array, BfSize i, BfConstPtr eltPtr) {
  BF_ERROR_BEGIN();

  BfSize n = array->size;

  if (n == array->capacity) {
    expand(array);
    HANDLE_ERROR();
  }

  bfMemMove(ELT_PTR(array, i), n - i, array->eltSize, ELT_PTR(array, i + 1));

  bfMemCopy(eltPtr, 1, array->eltSize, ELT_PTR(array, i));

  ++array->size;

  BF_ERROR_END() {
    BF_DIE();
  }
}

void bfArraySave(BfArray const *array, char const *path) {
  BF_ERROR_BEGIN();

  FILE *fp = fopen(path, "w");
  if (fp == NULL)
    RAISE_ERROR(BF_ERROR_FILE_ERROR);

  fwrite(array->data, array->eltSize, array->size, fp);

  BF_ERROR_END() {
    BF_DIE();
  }

  fclose(fp);
}

bool bfArrayIsEmpty(BfArray const *array) {
  return array->size == 0;
}

void bfArrayRemove(BfArray *array, BfSize i) {
  BF_ERROR_BEGIN();

  BfSize n = array->size;

  if (i >= n)
    RAISE_ERROR(BF_ERROR_OUT_OF_RANGE);

  bfMemMove(ELT_PTR(array, i + 1), n - i - 1, array->eltSize, ELT_PTR(array, i));

  --array->size;

  BF_ERROR_END() {
    BF_DIE();
  }
}
