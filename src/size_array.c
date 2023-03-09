#include <bf/size_array.h>

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <bf/error.h>
#include <bf/error_macros.h>

static BfSize const BF_ARRAY_DEFAULT_CAPACITY = 128;

static void invalidate(BfSizeArray *sizeArray) {
  sizeArray->data = NULL;
  sizeArray->size = BF_SIZE_BAD_VALUE;
  sizeArray->capacity = BF_SIZE_BAD_VALUE;
}

BfSizeArray *bfSizeArrayNew() {
  BEGIN_ERROR_HANDLING();

  BfSizeArray *sizeArray = malloc(sizeof(BfSizeArray));
  if (sizeArray == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  invalidate(sizeArray);

  END_ERROR_HANDLING() {
    sizeArray = NULL;
  }

  return sizeArray;
}

void bfSizeArrayInitWithDefaultCapacity(BfSizeArray *sizeArray) {
  BEGIN_ERROR_HANDLING();

  sizeArray->size = 0;
  sizeArray->capacity = BF_ARRAY_DEFAULT_CAPACITY;

  sizeArray->data = malloc(sizeArray->capacity*sizeof(BfSize));
  if (sizeArray->data == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING() {
    invalidate(sizeArray);
  }
}

void bfSizeArrayDeinit(BfSizeArray *sizeArray) {
  free(sizeArray->data);
  invalidate(sizeArray);
}

void bfSizeArrayDealloc(BfSizeArray **sizeArray) {
  free(*sizeArray);
  *sizeArray = NULL;
}

void bfSizeArrayDeinitAndDealloc(BfSizeArray **sizeArray) {
  bfSizeArrayDeinit(*sizeArray);
  bfSizeArrayDealloc(sizeArray);
}

void bfSizeArrayExpandCapacity(BfSizeArray *sizeArray, BfSize newCapacity) {
  BEGIN_ERROR_HANDLING();

  BfSize *data = NULL;

  if (newCapacity < sizeArray->capacity)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  data = malloc(newCapacity*sizeof(BfSize));
  if (data == NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  sizeArray->data = data;
  sizeArray->capacity = newCapacity;

  END_ERROR_HANDLING() {
    free(data);
  }
}

void bfSizeArrayAppend(BfSizeArray *sizeArray, BfSize elt) {
  BEGIN_ERROR_HANDLING();

  if (sizeArray->size == sizeArray->capacity) {
    bfSizeArrayExpandCapacity(sizeArray, 2*sizeArray->capacity);
    HANDLE_ERROR();
  }

  sizeArray->data[sizeArray->size++] = elt;

  END_ERROR_HANDLING() {
    assert(false);
  }
}

bool bfSizeArrayContains(BfSizeArray const *sizeArray, BfSize elt) {
  for (BfSize i = 0; i < sizeArray->size; ++i) {
    if (sizeArray->data[i] == elt)
      return true;
  }
  return false;
}

bool bfSizeArrayIsSorted(BfSizeArray const *sizeArray) {
  for (BfSize i = 1; i < sizeArray->size; ++i) {
    if (sizeArray->data[i - 1] > sizeArray->data[i])
      return false;
  }
  return true;
}

void bfSizeArrayInsertSorted(BfSizeArray *sizeArray, BfSize elt) {
  BEGIN_ERROR_HANDLING();

#if BF_DEBUG
  assert(bfSizeArrayIsSorted(sizeArray));

  BfSize prevSize = sizeArray->size;
#endif

  if (sizeArray->size == sizeArray->capacity) {
    bfSizeArrayExpandCapacity(sizeArray, 2*sizeArray->capacity);
    HANDLE_ERROR();
  }

  BfSize i = 0;
  for (; i < sizeArray->size; ++i)
    if (sizeArray->data[i] >= elt)
      break;

  memmove(sizeArray->data + i + 1, sizeArray->data + i,
          (sizeArray->size - i)*sizeof(BfSize));

  sizeArray->data[i] = elt;

  ++sizeArray->size;

  END_ERROR_HANDLING() {
    assert(false);
  }

#if BF_DEBUG
  assert(sizeArray->size == prevSize + 1);
#endif
}

BfSize bfSizeArrayFindFirst(BfSizeArray const *sizeArray, BfSize elt) {
  BfSize index = BF_SIZE_BAD_VALUE;
  for (BfSize i = 0; i < sizeArray->size; ++i) {
    if (sizeArray->data[i] == elt) {
      index = i;
      break;
    }
  }
  return index;
}

BfSize bfSizeArrayGet(BfSizeArray const *sizeArray, BfSize i) {
  BEGIN_ERROR_HANDLING();

  BfSize elt = BF_SIZE_BAD_VALUE;

  if (i >= sizeArray->size)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  elt = sizeArray->data[i];

  END_ERROR_HANDLING() {
    elt = BF_SIZE_BAD_VALUE;
  }

  return elt;
}

void bfSizeArrayCopyData(BfSizeArray const *sizeArray, BfSize *dst) {
  BEGIN_ERROR_HANDLING();

  if (dst == NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  memcpy(dst, sizeArray->data, sizeArray->size*sizeof(BfSize));

  END_ERROR_HANDLING() {}
}

BfSize bfSizeArrayGetSize(BfSizeArray const *sizeArray) {
  return sizeArray->size;
}
