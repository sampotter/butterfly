#include <bf/real_array.h>

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <bf/const.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/rand.h>
#include <bf/vec_real.h>

static void invalidate(BfRealArray *realArray) {
  realArray->data = NULL;
  realArray->size = BF_SIZE_BAD_VALUE;
  realArray->capacity = BF_SIZE_BAD_VALUE;
}

BfRealArray *bfRealArrayNew() {
  BEGIN_ERROR_HANDLING();

  BfRealArray *realArray = malloc(sizeof(BfRealArray));
  if (realArray == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  invalidate(realArray);

  END_ERROR_HANDLING() {
    realArray = NULL;
  }

  return realArray;
}

BfRealArray *bfRealArrayNewWithDefaultCapacity() {
  BEGIN_ERROR_HANDLING();

  BfRealArray *realArray = bfRealArrayNew();
  HANDLE_ERROR();

  bfRealArrayInitWithDefaultCapacity(realArray);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {
    assert(false);
  }

  return realArray;
}

void bfRealArrayInitWithDefaultCapacity(BfRealArray *realArray) {
  BEGIN_ERROR_HANDLING();

  realArray->size = 0;
  realArray->capacity = BF_ARRAY_DEFAULT_CAPACITY;

  realArray->data = malloc(realArray->capacity*sizeof(BfReal));
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
  free(realArray->data);
  invalidate(realArray);
}

void bfRealArrayDealloc(BfRealArray **realArray) {
  free(*realArray);
  *realArray = NULL;
}

void bfRealArrayDeinitAndDealloc(BfRealArray **realArray) {
  bfRealArrayDeinit(*realArray);
  bfRealArrayDealloc(realArray);
}

void bfRealArrayExpandCapacity(BfRealArray *realArray, BfSize newCapacity) {
  BEGIN_ERROR_HANDLING();

  BfReal *data = NULL, *oldData = NULL;

  if (newCapacity < realArray->capacity)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  data = malloc(newCapacity*sizeof(BfReal));
  if (data == NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

#if BF_DEBUG
  for (BfSize i = 0; i < newCapacity; ++i)
    data[i] = BF_NAN;
#endif

  /* Copy over old values */
  memcpy(data, realArray->data, realArray->size*sizeof(BfReal));

  oldData = realArray->data;

  realArray->data = data;
  realArray->capacity = newCapacity;

  END_ERROR_HANDLING() {
    assert(false);
  }

  free(oldData);
}

void bfRealArrayAppend(BfRealArray *realArray, BfReal elt) {
  BEGIN_ERROR_HANDLING();

  if (realArray->size == realArray->capacity) {
    bfRealArrayExpandCapacity(realArray, 2*realArray->capacity);
    HANDLE_ERROR();
  }

  realArray->data[realArray->size++] = elt;

  END_ERROR_HANDLING() {
    assert(false);
  }
}

BfVec *bfRealArrayGetVecView(BfRealArray *realArray) {
  BEGIN_ERROR_HANDLING();

  BfVecReal *vecReal = bfVecRealNew();
  HANDLE_ERROR();

  bfVecRealInitView(vecReal, realArray->size, 1, realArray->data);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {
    assert(false);
  }

  return bfVecRealToVec(vecReal);
}
