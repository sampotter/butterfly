#include <bf/size_array.h>

#include <bf/assert.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mem.h>
#include <bf/rand.h>
#include <bf/util.h>

static void invalidate(BfSizeArray *sizeArray) {
  sizeArray->data = NULL;
  sizeArray->size = BF_SIZE_BAD_VALUE;
  sizeArray->capacity = BF_SIZE_BAD_VALUE;
}

BfSizeArray *bfSizeArrayNew() {
  BF_ERROR_BEGIN();

  BfSizeArray *sizeArray = bfMemAlloc(1, sizeof(BfSizeArray));
  if (sizeArray == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  invalidate(sizeArray);

  BF_ERROR_END() {
    sizeArray = NULL;
  }

  return sizeArray;
}

BfSizeArray *bfSizeArrayNewIota(BfSize n) {
  BF_ERROR_BEGIN();

  BfSizeArray *sizeArray = bfSizeArrayNewWithDefaultCapacity();
  HANDLE_ERROR();

  for (BfSize i = 0; i < n; ++i)
    bfSizeArrayAppend(sizeArray, i);

  BF_ERROR_END() {
    BF_DIE();
  }

  return sizeArray;
}

BfSizeArray *bfSizeArrayNewWithCapacity(BfSize capacity) {
  BF_ERROR_BEGIN();

  BfSizeArray *sizeArray = bfSizeArrayNew();
  HANDLE_ERROR();

  bfSizeArrayInitWithCapacity(sizeArray, capacity);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  return sizeArray;
}

BfSizeArray *bfSizeArrayNewWithDefaultCapacity() {
  return bfSizeArrayNewWithCapacity(BF_ARRAY_DEFAULT_CAPACITY);
}

void bfSizeArrayInitWithCapacity(BfSizeArray *sizeArray, BfSize capacity) {
  BF_ERROR_BEGIN();

  sizeArray->size = 0;
  sizeArray->capacity = capacity;

  sizeArray->data = bfMemAlloc(sizeArray->capacity, sizeof(BfSize));
  if (sizeArray->data == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

#if BF_DEBUG
  for (BfSize i = 0; i < sizeArray->capacity; ++i)
    sizeArray->data[i] = BF_SIZE_BAD_VALUE;
#endif

  BF_ERROR_END() {
    invalidate(sizeArray);
  }
}

void bfSizeArrayInitWithDefaultCapacity(BfSizeArray *sizeArray) {
  bfSizeArrayInitWithCapacity(sizeArray, BF_ARRAY_DEFAULT_CAPACITY);
}

void bfSizeArrayDeinit(BfSizeArray *sizeArray) {
  bfMemFree(sizeArray->data);
  invalidate(sizeArray);
}

void bfSizeArrayDealloc(BfSizeArray **sizeArray) {
  bfMemFree(*sizeArray);
  *sizeArray = NULL;
}

void bfSizeArrayDeinitAndDealloc(BfSizeArray **sizeArray) {
  bfSizeArrayDeinit(*sizeArray);
  bfSizeArrayDealloc(sizeArray);
}

void bfSizeArrayExpandCapacity(BfSizeArray *sizeArray, BfSize newCapacity) {
  BF_ERROR_BEGIN();

  BfSize *data = NULL, *oldData = NULL;

  if (newCapacity < sizeArray->capacity)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  data = bfMemAlloc(newCapacity, sizeof(BfSize));
  if (data == NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

#if BF_DEBUG
  for (BfSize i = 0; i < newCapacity; ++i)
    data[i] = BF_SIZE_BAD_VALUE;
#endif

  /* Copy over old values */
  bfMemCopy(sizeArray->data, sizeArray->size, sizeof(BfSize), data);

  oldData = sizeArray->data;

  sizeArray->data = data;
  sizeArray->capacity = newCapacity;

  BF_ERROR_END() {
    BF_DIE();
  }

  bfMemFree(oldData);
}

void bfSizeArrayAppend(BfSizeArray *sizeArray, BfSize elt) {
  BF_ERROR_BEGIN();

  if (elt == BF_SIZE_BAD_VALUE)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (sizeArray->size == sizeArray->capacity) {
    bfSizeArrayExpandCapacity(sizeArray, 2*sizeArray->capacity);
    HANDLE_ERROR();
  }

  sizeArray->data[sizeArray->size++] = elt;

  BF_ERROR_END() {
    BF_DIE();
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

static int sizeCompar(const void *elt1, const void *elt2, void *_) {
  ComparatorAndAuxPtr *cmpAndAux = (ComparatorAndAuxPtr *)_;
  BfSizeArrayComparator cmp = cmpAndAux->cmp;
  void *aux = cmpAndAux->aux;
  return cmp(*(BfSize const *)elt1, *(BfSize const *)elt2, aux);
}

void bfSizeArraySort(BfSizeArray *sizeArray, BfSizeArrayComparator cmp, void *aux) {
  ComparatorAndAuxPtr cmpAndAux = {.cmp = cmp, .aux = aux};
  bfSort(sizeArray->data, sizeArray->size, sizeof(BfSize), (BfCompar)sizeCompar, &cmpAndAux);
}

void bfSizeArrayInsertSorted(BfSizeArray *sizeArray, BfSize elt) {
  BF_ERROR_BEGIN();

#if BF_DEBUG
  BF_ASSERT(bfSizeArrayIsSorted(sizeArray));

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

  bfMemMove(sizeArray->data + i, sizeArray->size - i, sizeof(BfSize), sizeArray->data + i + 1);

  sizeArray->data[i] = elt;

  ++sizeArray->size;

  BF_ERROR_END() {
    BF_DIE();
  }

#if BF_DEBUG
  BF_ASSERT(sizeArray->size == prevSize + 1);
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
  BF_ERROR_BEGIN();

  BfSize elt = BF_SIZE_BAD_VALUE;

  if (i >= sizeArray->size)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  elt = sizeArray->data[i];

  BF_ERROR_END() {
    elt = BF_SIZE_BAD_VALUE;
  }

  return elt;
}

BfSize bfSizeArrayGetFirst(BfSizeArray const *sizeArray) {
  BF_ERROR_BEGIN();

  BfSize elt = BF_SIZE_BAD_VALUE;

  if (sizeArray->size == 0)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  elt = sizeArray->data[0];

  BF_ERROR_END() {
    BF_DIE();
  }

  return elt;
}

BfSize bfSizeArrayGetLast(BfSizeArray const *sizeArray) {
  BF_ERROR_BEGIN();

  BfSize elt = BF_SIZE_BAD_VALUE;

  if (sizeArray->size == 0)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  elt = sizeArray->data[sizeArray->size - 1];

  BF_ERROR_END() {
    BF_DIE();
  }

  return elt;
}

BfSize bfSizeArrayGetRand(BfSizeArray const *sizeArray) {
  return sizeArray->size == 0 ?
    BF_SIZE_BAD_VALUE :
    sizeArray->data[bfSizeUniform1(0, sizeArray->size)];
}

void bfSizeArrayCopyData(BfSizeArray const *sizeArray, BfSize *dst) {
  BF_ERROR_BEGIN();

  if (dst == NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  bfMemCopy(sizeArray->data, sizeArray->size, sizeof(BfSize), dst);

  BF_ERROR_END() {}
}

BfSize bfSizeArrayGetSize(BfSizeArray const *sizeArray) {
  return sizeArray->size;
}

void bfSizeArrayDelete(BfSizeArray *sizeArray, BfSize i) {
  BF_ERROR_BEGIN();

  BfSize n = sizeArray->size;

  if (i >= n)
    RAISE_ERROR(BF_ERROR_OUT_OF_RANGE);

  BfSize *ptr = sizeArray->data;

  bfMemMove(ptr + i + 1, n - i - 1, sizeof(BfSize), ptr + i);

  --sizeArray->size;

  BF_ERROR_END() {
    BF_DIE();
  }
}

void bfSizeArrayDeleteFirst(BfSizeArray *sizeArray, BfSize elt) {
  BfSize i = bfSizeArrayFindFirst(sizeArray, elt);

  if (i == BF_SIZE_BAD_VALUE)
    return;

  BfSize *dst = sizeArray->data + i;
  BfSize *src = dst + 1;
  bfMemMove(src, sizeArray->size - i, sizeof(BfSize), dst);

  --sizeArray->size;
}

void bfSizeArraySave(BfSizeArray const *sizeArray, char const *path) {
  BF_ERROR_BEGIN();

  FILE *fp = fopen(path, "w");
  if (fp == NULL)
    RAISE_ERROR(BF_ERROR_FILE_ERROR);

  fwrite(sizeArray->data, sizeof(BfSize), sizeArray->size, fp);

  fclose(fp);

  BF_ERROR_END() {
    BF_DIE();
  }
}

BfSize *bfSizeArrayGetDataPtr(BfSizeArray const *sizeArray) {
  return sizeArray->data;
}
