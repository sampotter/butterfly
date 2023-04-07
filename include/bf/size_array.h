#pragma once

#include "def.h"
#include "types.h"

struct BfSizeArray {
  BfSize *data;
  BfSize size;
  BfSize capacity;
};

typedef int (*BfSizeArrayComparator)(BfSize, BfSize, void *);

BfSizeArray *bfSizeArrayNew();
BfSizeArray *bfSizeArrayNewIota(BfSize n);
BfSizeArray *bfSizeArrayNewWithDefaultCapacity();
void bfSizeArrayInitWithDefaultCapacity(BfSizeArray *sizeArray);
void bfSizeArrayDeinit(BfSizeArray *sizeArray);
void bfSizeArrayDealloc(BfSizeArray **sizeArray);
void bfSizeArrayDeinitAndDealloc(BfSizeArray **sizeArray);
void bfSizeArrayExpandCapacity(BfSizeArray *sizeArray, BfSize newCapacity);
void bfSizeArrayAppend(BfSizeArray *sizeArray, BfSize elt);
bool bfSizeArrayContains(BfSizeArray const *sizeArray, BfSize elt);
bool bfSizeArrayIsEmpty(BfSizeArray const *sizeArray);
bool bfSizeArrayIsSorted(BfSizeArray const *sizeArray);
void bfSizeArraySort(BfSizeArray *sizeArray, BfSizeArrayComparator cmp, void *aux);
void bfSizeArrayInsertSorted(BfSizeArray *sizeArray, BfSize elt);
BfSize bfSizeArrayFindFirst(BfSizeArray const *sizeArray, BfSize elt);
BfSize bfSizeArrayGet(BfSizeArray const *sizeArray, BfSize i);
BfSize bfSizeArrayGetLast(BfSizeArray const *sizeArray);
BfSize bfSizeArrayGetRand(BfSizeArray const *sizeArray);
void bfSizeArrayCopyData(BfSizeArray const *sizeArray, BfSize *dst);
BfSize bfSizeArrayGetSize(BfSizeArray const *sizeArray);
void bfSizeArrayDeleteFirst(BfSizeArray *sizeArray, BfSize elt);
