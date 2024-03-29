#include <bf/real_array.h>

#include <bf/assert.h>
#include <bf/const.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mem.h>
#include <bf/rand.h>
#include <bf/vec_real.h>

#include "macros.h"

static void invalidate(BfRealArray *realArray) {
  realArray->data = NULL;
  realArray->size = BF_SIZE_BAD_VALUE;
  realArray->capacity = BF_SIZE_BAD_VALUE;
  realArray->isView = false;
}

BfRealArray *bfRealArrayCopy(BfRealArray const *realArray) {
  BF_ERROR_BEGIN();

  BfRealArray *realArrayCopy = bfMemAlloc(1, sizeof(BfRealArray));
  HANDLE_ERROR();

  invalidate(realArrayCopy);

  bfRealArrayInitCopy(realArrayCopy, realArray);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  return realArrayCopy;
}

BfRealArray *bfRealArrayNew() {
  BF_ERROR_BEGIN();

  BfRealArray *realArray = bfMemAlloc(1, sizeof(BfRealArray));
  if (realArray == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  invalidate(realArray);

  BF_ERROR_END() {
    realArray = NULL;
  }

  return realArray;
}

BfRealArray *bfRealArrayNewFromVecReal(BfVecReal const *vecReal, BfPolicy policy) {
  BF_ERROR_BEGIN();

  BfRealArray *realArray = bfRealArrayNew();
  HANDLE_ERROR();

  if (policy == BF_POLICY_VIEW) {
    realArray->size = vecReal->super.size;
    realArray->capacity = BF_SIZE_BAD_VALUE;
    realArray->isView = true;
    realArray->data = vecReal->data;
  }

  else if (policy == BF_POLICY_COPY) {
    realArray->size = vecReal->super.size;
    realArray->capacity = vecReal->super.size;
    realArray->isView = false;

    realArray->data = bfMemAlloc(realArray->size, sizeof(BfReal));
    HANDLE_ERROR();

    BfReal const *readPtr = vecReal->data;
    BfReal *writePtr = realArray->data;
    for (BfSize i = 0; i < realArray->size; ++i) {
      *writePtr = *readPtr;
      readPtr += vecReal->stride;
      ++writePtr;
    }
  }

  else RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);

  BF_ERROR_END() {
    BF_DIE();
  }

  return realArray;
}

BfRealArray *bfRealArrayNewWithCapacity(BfSize capacity) {
  BF_ERROR_BEGIN();

  BfRealArray *realArray = bfRealArrayNew();
  HANDLE_ERROR();

  bfRealArrayInitWithCapacity(realArray, capacity);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  return realArray;
}

BfRealArray *bfRealArrayNewWithDefaultCapacity(void) {
  BF_ERROR_BEGIN();

  BfRealArray *realArray = bfRealArrayNew();
  HANDLE_ERROR();

  bfRealArrayInitWithDefaultCapacity(realArray);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  return realArray;
}

BfRealArray *bfRealArrayNewWithValue(BfSize size, BfReal value) {
  BF_ERROR_BEGIN();

  BfRealArray *realArray = bfRealArrayNew();
  HANDLE_ERROR();

  bfRealArrayInitWithValue(realArray, size, value);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  return realArray;
}

BfRealArray *bfRealArrayNewFromFile(char const *path) {
  BF_ERROR_BEGIN();

  FILE *fp = fopen(path, "r");
  if (fp == NULL)
    RAISE_ERROR(BF_ERROR_FILE_ERROR);

  fseek(fp, 0, SEEK_END);
  BfSize numFileBytes = ftell(fp);
  fseek(fp, 0, SEEK_SET);

  if (numFileBytes % sizeof(BfReal) != 0)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  BfSize capacity = numFileBytes/sizeof(BfReal);

  BfRealArray *realArray = bfRealArrayNewWithCapacity(capacity);
  HANDLE_ERROR();

  realArray->size = capacity;
  fread(realArray->data, sizeof(BfReal), realArray->size, fp);

  BF_ERROR_END() {
    BF_DIE();
  }

  fclose(fp);

  return realArray;
}

void bfRealArrayInitWithCapacity(BfRealArray *realArray, BfSize capacity) {
  BF_ERROR_BEGIN();

  realArray->size = 0;
  realArray->capacity = capacity;
  realArray->isView = false;

  realArray->data = bfMemAlloc(realArray->capacity, sizeof(BfReal));
  if (realArray->data == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

#if BF_DEBUG
  for (BfSize i = 0; i < realArray->capacity; ++i)
    realArray->data[i] = BF_NAN;
#endif

  BF_ERROR_END() {
    invalidate(realArray);
  }
}

void bfRealArrayInitWithDefaultCapacity(BfRealArray *realArray) {
  bfRealArrayInitWithCapacity(realArray, BF_ARRAY_DEFAULT_CAPACITY);
}

void bfRealArrayInitWithValue(BfRealArray *realArray, BfSize size, BfReal value) {
  BF_ERROR_BEGIN();

  realArray->size = size;
  realArray->capacity = size;
  realArray->isView = false;

  realArray->data = bfMemAlloc(realArray->capacity, sizeof(BfReal));
  HANDLE_ERROR();

  for (BfSize i = 0; i < realArray->size; ++i) {
    realArray->data[i] = value;
  }

  BF_ERROR_END() {
    BF_DIE();
  }
}

void bfRealArrayInitCopy(BfRealArray *realArray, BfRealArray const *otherRealArray) {
  BF_ERROR_BEGIN();

  realArray->size = otherRealArray->size;
  realArray->capacity = otherRealArray->size;
  realArray->isView = false;

  realArray->data = bfMemAlloc(realArray->size, sizeof(BfReal));
  HANDLE_ERROR();

  bfMemCopy(otherRealArray->data, otherRealArray->size, sizeof(BfSize), realArray->data);

  BF_ERROR_END() {
    BF_DIE();
  }
}

void bfRealArrayDeinit(BfRealArray *realArray) {
  if (!realArray->isView)
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

  if (realArray->isView)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

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

  BF_ERROR_END() {
    BF_DIE();
  }

  bfMemFree(oldData);
}

void bfRealArrayShrinkCapacityToSize(BfRealArray *realArray) {
  BF_ERROR_BEGIN();

  BF_ASSERT(realArray->capacity >= realArray->size);

  bfMemRealloc(realArray->data, realArray->capacity, sizeof(BfReal));
  HANDLE_ERROR();

  realArray->capacity = realArray->size;

  BF_ERROR_END() {
    BF_DIE();
  }
}

void bfRealArrayAppend(BfRealArray *realArray, BfReal elt) {
  BF_ERROR_BEGIN();

  if (realArray->isView)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  if (realArray->size == realArray->capacity) {
    bfRealArrayExpandCapacity(realArray, 2*realArray->capacity);
    HANDLE_ERROR();
  }

  realArray->data[realArray->size++] = elt;

  BF_ERROR_END() {
    BF_DIE();
  }
}

void bfRealArrayExtend(BfRealArray *realArray, BfRealArray const *otherRealArray) {
  BF_ERROR_BEGIN();

  for (BfSize i = 0; i < bfRealArrayGetSize(otherRealArray); ++i) {
    bfRealArrayAppend(realArray, bfRealArrayGetValue(otherRealArray, i));
    HANDLE_ERROR();
  }

  BF_ERROR_END() {
    BF_DIE();
  }
}

BfVec *bfRealArrayGetVecView(BfRealArray *realArray) {
  BF_ERROR_BEGIN();

  BfVecReal *vecReal = bfVecRealNew();
  HANDLE_ERROR();

  bfVecRealInitView(vecReal, realArray->size, 1, realArray->data);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
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

  BF_ERROR_END() {
    BF_DIE();
  }

  return bfVecRealToVec(vecReal);
}

BfReal bfRealArrayGetValue(BfRealArray const *realArray, BfSize i) {
  BF_ERROR_BEGIN();

  BfReal value = BF_NAN;

  if (i == BF_SIZE_BAD_VALUE)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (i >= realArray->size)
    RAISE_ERROR(BF_ERROR_OUT_OF_RANGE);

  value = realArray->data[i];

  BF_ERROR_END() {
    BF_DIE();
  }

  return value;
}

void bfRealArraySetValue(BfRealArray *realArray, BfSize i, BfReal value) {
  BF_ERROR_BEGIN();

  if (i == BF_SIZE_BAD_VALUE)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (i >= realArray->size)
    RAISE_ERROR(BF_ERROR_OUT_OF_RANGE);

  realArray->data[i] = value;

  BF_ERROR_END() {
    BF_DIE();
  }
}

void bfRealArrayGetValues(BfRealArray const *realArray, BfSize n, BfSize const *inds, BfReal *values) {
  BF_ERROR_BEGIN();

  for (BfSize i = 0; i < n; ++i)
    if (inds[i] >= realArray->size)
      RAISE_ERROR(BF_ERROR_OUT_OF_RANGE);

  bfMemCopy(realArray->data, n, sizeof(BfReal), values);

  BF_ERROR_END() {
    BF_DIE();
  }
}

void bfRealArrayInsert(BfRealArray *realArray, BfSize i, BfReal value) {
  BF_ERROR_BEGIN();

  if (i > realArray->size)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (realArray->isView)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  if (realArray->size == realArray->capacity) {
    bfRealArrayExpandCapacity(realArray, 2*realArray->capacity);
    HANDLE_ERROR();
  }

  BfSize n = realArray->size;

  BfReal *newData = bfMemRealloc(realArray->data, n + 1, sizeof(BfReal));
  HANDLE_ERROR();

  bfMemMove(newData + i, n - i, sizeof(BfSize), newData + i + 1);

  newData[i] = value;

  realArray->data = newData;

  ++realArray->size;

  BF_ERROR_END() {
    BF_DIE();
  }
}

BfSize bfRealArrayGetSize(BfRealArray const *realArray) {
  return realArray->size;
}

bool bfRealArrayIsEmpty(BfRealArray const *realArray) {
  return realArray->size == 0;
}

void bfRealArraySave(BfRealArray const *realArray, char const *path) {
  BF_ERROR_BEGIN();

  FILE *fp = fopen(path, "w");
  if (fp == NULL)
    RAISE_ERROR(BF_ERROR_FILE_ERROR);

  fwrite(realArray->data, sizeof(BfReal), realArray->size, fp);

  BF_ERROR_END() {
    BF_DIE();
  }

  fclose(fp);
}

void bfRealArrayNegate(BfRealArray *realArray) {
  for (BfSize i = 0; i < realArray->size; ++i)
    realArray->data[i] *= -1;
}

BfPerm *bfRealArrayArgsort(BfRealArray const *realArray) {
  BF_ERROR_BEGIN();

  BfSize const n = bfRealArrayGetSize(realArray);

  BfPerm *perm = bfPermNewIdentity(n);
  HANDLE_ERROR();

  /* TODO: just using selection sort here for now... should upgrade to
   * something better later */
  if (n > 0) {
    for (BfSize i = 0; i < n - 1; ++i) {
      BfSize k = i;
      BfReal value_k = bfRealArrayGetValue(realArray, perm->index[k]);
      for (BfSize j = i + 1; j < n; ++j) {
        BfReal const value_j = bfRealArrayGetValue(realArray, perm->index[j]);
        if (value_j < value_k) {
          k = j;
          value_k = value_j;
        }
      }
      BF_SWAP(perm->index[i], perm->index[k]);
    }
  }

  bfPermReverse(perm);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  return perm;
}

void bfRealArrayPermute(BfRealArray *realArray, BfPerm const *perm) {
  BF_ERROR_BEGIN();

  BfRealArray *realArrayCopy = bfRealArrayCopy(realArray);
  HANDLE_ERROR();

  for (BfSize i = 0; i < realArray->size; ++i) {
    BfReal const value = bfRealArrayGetValue(realArrayCopy, i);
    bfRealArraySetValue(realArray, perm->index[i], value);
  }

  BF_ERROR_END() {
    BF_DIE();
  }

  bfRealArrayDeinitAndDealloc(&realArrayCopy);
}

BfReal *bfRealArrayGetDataPtr(BfRealArray *realArray) {
  return realArray->data;
}

BfReal *bfRealArrayStealPtr(BfRealArray *realArray) {
  BF_ERROR_BEGIN();

  if (realArray->isView)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  realArray->isView = true;

  BF_ERROR_END() {
    BF_DIE();
  }

  return realArray->data;
}
