#include <bf/real_array.h>

#include <bf/assert.h>
#include <bf/const.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mem.h>
#include <bf/rand.h>
#include <bf/vec_real.h>

static void invalidate(BfRealArray *realArray) {
  realArray->data = NULL;
  realArray->size = BF_SIZE_BAD_VALUE;
  realArray->capacity = BF_SIZE_BAD_VALUE;
}

BfRealArray *bfRealArrayNew() {
  BF_ERROR_BEGIN();

  BfRealArray *realArray = bfMemAlloc(1, sizeof(BfRealArray));
  if (realArray == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  invalidate(realArray);

  END_ERROR_HANDLING() {
    realArray = NULL;
  }

  return realArray;
}

BfRealArray *bfRealArrayNewWithDefaultCapacity() {
  BF_ERROR_BEGIN();

  BfRealArray *realArray = bfRealArrayNew();
  HANDLE_ERROR();

  bfRealArrayInitWithDefaultCapacity(realArray);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {
    BF_ASSERT(false);
  }

  return realArray;
}

void bfRealArrayInitWithDefaultCapacity(BfRealArray *realArray) {
  BF_ERROR_BEGIN();

  realArray->size = 0;
  realArray->capacity = BF_ARRAY_DEFAULT_CAPACITY;

  realArray->data = bfMemAlloc(realArray->capacity, sizeof(BfReal));
  if (realArray->data == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

#if BF_DEBUG
  for (BfSize i = 0; i < realArray->capacity; ++i)
    realArray->data[i] = BF_NAN;
#endif

  END_ERROR_HANDLING() {
    invalidate(realArray);
  }
}

void bfRealArrayDeinit(BfRealArray *realArray) {
  bfMemFree(realArray->data);
  invalidate(realArray);
}

void bfRealArrayDealloc(BfRealArray **realArray) {
  bfMemFree(*realArray);
  *realArray = NULL;
}

void bfRealArrayDeinitAndDealloc(BfRealArray **realArray) {
  bfRealArrayDeinit(*realArray);
  bfRealArrayDealloc(realArray);
}

void bfRealArrayExpandCapacity(BfRealArray *realArray, BfSize newCapacity) {
  BF_ERROR_BEGIN();

  BfReal *data = NULL, *oldData = NULL;

  if (newCapacity < realArray->capacity)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  data = bfMemAlloc(newCapacity, sizeof(BfReal));
  if (data == NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

#if BF_DEBUG
  for (BfSize i = 0; i < newCapacity; ++i)
    data[i] = BF_NAN;
#endif

  /* Copy over old values */
  bfMemCopy(realArray->data, realArray->size, sizeof(BfReal), data);

  oldData = realArray->data;

  realArray->data = data;
  realArray->capacity = newCapacity;

  END_ERROR_HANDLING() {
    BF_ASSERT(false);
  }

  bfMemFree(oldData);
}

void bfRealArrayAppend(BfRealArray *realArray, BfReal elt) {
  BF_ERROR_BEGIN();

  if (realArray->size == realArray->capacity) {
    bfRealArrayExpandCapacity(realArray, 2*realArray->capacity);
    HANDLE_ERROR();
  }

  realArray->data[realArray->size++] = elt;

  END_ERROR_HANDLING() {
    BF_ASSERT(false);
  }
}

BfVec *bfRealArrayGetVecView(BfRealArray *realArray) {
  BF_ERROR_BEGIN();

  BfVecReal *vecReal = bfVecRealNew();
  HANDLE_ERROR();

  bfVecRealInitView(vecReal, realArray->size, 1, realArray->data);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {
    BF_ASSERT(false);
  }

  return bfVecRealToVec(vecReal);
}

BfVec *bfRealArrayGetSubvecView(BfRealArray *realArray, BfSize i0, BfSize i1) {
  BF_ERROR_BEGIN();

  if (i0 > i1)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (i1 > realArray->size)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize size = i1 - i0;

  BfVecReal *vecReal = bfVecRealNew();
  HANDLE_ERROR();

  bfVecRealInitView(vecReal, size, 1, realArray->data + i0);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {
    BF_ASSERT(false);
  }

  return bfVecRealToVec(vecReal);
}
