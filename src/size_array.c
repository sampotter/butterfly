#include <bf/size_array.h>

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/rand.h>

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

BfSizeArray *bfSizeArrayNewIota(BfSize n) {
  BEGIN_ERROR_HANDLING();

  BfSizeArray *sizeArray = bfSizeArrayNewWithDefaultCapacity();
  HANDLE_ERROR();

  for (BfSize i = 0; i < n; ++i)
    bfSizeArrayAppend(sizeArray, i);

  END_ERROR_HANDLING() {
    assert(false);
  }

  return sizeArray;
}

BfSizeArray *bfSizeArrayNewWithDefaultCapacity() {
  BEGIN_ERROR_HANDLING();

  BfSizeArray *sizeArray = bfSizeArrayNew();
  HANDLE_ERROR();

  bfSizeArrayInitWithDefaultCapacity(sizeArray);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {
    assert(false);
  }

  return sizeArray;
}

void bfSizeArrayInitWithDefaultCapacity(BfSizeArray *sizeArray) {
  BEGIN_ERROR_HANDLING();

  sizeArray->size = 0;
  sizeArray->capacity = 128;

  sizeArray->data = malloc(sizeArray->capacity*sizeof(BfSize));
  if (sizeArray->data == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

#if BF_DEBUG
  for (BfSize i = 0; i < sizeArray->capacity; ++i)
    sizeArray->data[i] = BF_SIZE_BAD_VALUE;
#endif

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

  BfSize *data = NULL, *oldData = NULL;

  if (newCapacity < sizeArray->capacity)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  data = malloc(newCapacity*sizeof(BfSize));
  if (data == NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

#if BF_DEBUG
  for (BfSize i = 0; i < newCapacity; ++i)
    data[i] = BF_SIZE_BAD_VALUE;
#endif

  /* Copy over old values */
  memcpy(data, sizeArray->data, sizeArray->size*sizeof(BfSize));

  oldData = sizeArray->data;

  sizeArray->data = data;
  sizeArray->capacity = newCapacity;

  END_ERROR_HANDLING() {
    assert(false);
  }

  free(oldData);
}

void bfSizeArrayAppend(BfSizeArray *sizeArray, BfSize elt) {
  BEGIN_ERROR_HANDLING();

  if (elt == BF_SIZE_BAD_VALUE)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

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

bool bfSizeArrayIsEmpty(BfSizeArray const *sizeArray) {
  return sizeArray->size == 0;
}

bool bfSizeArrayIsSorted(BfSizeArray const *sizeArray) {
  for (BfSize i = 1; i < sizeArray->size; ++i) {
    if (sizeArray->data[i - 1] > sizeArray->data[i])
      return false;
  }
  return true;
}

typedef struct {
  BfSizeArrayComparator cmp;
  void *aux;
} ComparatorAndAuxPtr;

static int qsortComparator(const void *elt1, const void *elt2, void *_) {
  ComparatorAndAuxPtr *cmpAndAux = (ComparatorAndAuxPtr *)_;
  BfSizeArrayComparator cmp = cmpAndAux->cmp;
  void *aux = cmpAndAux->aux;
  return cmp(*(BfSize const *)elt1, *(BfSize const *)elt2, aux);
}

void bfSizeArraySort(BfSizeArray *sizeArray, BfSizeArrayComparator cmp, void *aux) {
  ComparatorAndAuxPtr cmpAndAux = {.cmp = cmp, .aux = aux};
  qsort_r(sizeArray->data, sizeArray->size, sizeof(BfSize), qsortComparator, &cmpAndAux);
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

BfSize bfSizeArrayGetFirst(BfSizeArray const *sizeArray) {
  BEGIN_ERROR_HANDLING();

  BfSize elt = BF_SIZE_BAD_VALUE;

  if (sizeArray->size == 0)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  elt = sizeArray->data[0];

  END_ERROR_HANDLING() {
    assert(false);
  }

  return elt;
}

BfSize bfSizeArrayGetLast(BfSizeArray const *sizeArray) {
  BEGIN_ERROR_HANDLING();

  BfSize elt = BF_SIZE_BAD_VALUE;

  if (sizeArray->size == 0)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  elt = sizeArray->data[sizeArray->size - 1];

  END_ERROR_HANDLING() {
    assert(false);
  }

  return elt;
}

BfSize bfSizeArrayGetRand(BfSizeArray const *sizeArray) {
  return sizeArray->size == 0 ?
    BF_SIZE_BAD_VALUE :
    sizeArray->data[bfSizeUniform1(0, sizeArray->size)];
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

void bfSizeArrayDeleteFirst(BfSizeArray *sizeArray, BfSize elt) {
  BfSize i = bfSizeArrayFindFirst(sizeArray, elt);

  if (i == BF_SIZE_BAD_VALUE)
    return;

  BfSize *dst = sizeArray->data + i;
  BfSize *src = dst + 1;
  BfSize size = (sizeArray->size - i)*sizeof(BfSize);

  memmove(dst, src, size);

  --sizeArray->size;
}
