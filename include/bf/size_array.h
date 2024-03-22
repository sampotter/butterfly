#pragma once

#include "def.h"
#include "types.h"

struct BfSizeArray {
  BfSize *data;
  BfSize size;
  BfSize capacity;
  bool isView;
};

typedef int (*BfSizeArrayComparator)(BfSize, BfSize, void *);

BfSizeArray *bfSizeArrayNew(void);
BfSizeArray *bfSizeArrayNewIota(BfSize n);
BfSizeArray *bfSizeArrayNewWithCapacity(BfSize capacity);
BfSizeArray *bfSizeArrayNewWithDefaultCapacity(void);
BfSizeArray *bfSizeArrayNewView(BfSize n, BfSize *inds);
BfSizeArray *bfSizeArrayNewFromFile(char const *path);
void bfSizeArrayInitWithCapacity(BfSizeArray *sizeArray, BfSize capacity);
void bfSizeArrayInitWithDefaultCapacity(BfSizeArray *sizeArray);
void bfSizeArrayInitView(BfSizeArray *sizeArray, BfSize n, BfSize *inds);
void bfSizeArrayDeinit(BfSizeArray *sizeArray);
void bfSizeArrayDealloc(BfSizeArray **sizeArray);
void bfSizeArrayDeinitAndDealloc(BfSizeArray **sizeArray);
void bfSizeArrayExpandCapacity(BfSizeArray *sizeArray, BfSize newCapacity);
void bfSizeArrayShrinkCapacityToSize(BfSizeArray *realArray);
void bfSizeArrayAppend(BfSizeArray *sizeArray, BfSize elt);
void bfSizeArrayExtend(BfSizeArray *sizeArray, BfSizeArray const *otherSizeArray);
bool bfSizeArrayContains(BfSizeArray const *sizeArray, BfSize elt);
bool bfSizeArrayIsEmpty(BfSizeArray const *sizeArray);
bool bfSizeArrayIsSorted(BfSizeArray const *sizeArray);
void bfSizeArraySort(BfSizeArray *sizeArray, BfSizeArrayComparator cmp, void *aux);
void bfSizeArrayInsertSorted(BfSizeArray *sizeArray, BfSize elt);
BfSize bfSizeArrayFindFirst(BfSizeArray const *sizeArray, BfSize elt);
BfSize bfSizeArrayGet(BfSizeArray const *sizeArray, BfSize i);
BfSize bfSizeArrayGetFirst(BfSizeArray const *sizeArray);
BfSize bfSizeArrayGetLast(BfSizeArray const *sizeArray);
BfSize bfSizeArrayGetRand(BfSizeArray const *sizeArray);
void bfSizeArrayCopyData(BfSizeArray const *sizeArray, BfSize *dst);
BfSize bfSizeArrayGetSize(BfSizeArray const *sizeArray);
void bfSizeArrayDelete(BfSizeArray *sizeArray, BfSize i);
void bfSizeArrayDeleteFirst(BfSizeArray *sizeArray, BfSize elt);
void bfSizeArraySave(BfSizeArray const *sizeArray, char const *path);
BfSize *bfSizeArrayGetDataPtr(BfSizeArray const *sizeArray);
BfSize *bfSizeArrayStealPtr(BfSizeArray *sizeArray);
BfSize bfSizeArrayGetMaximum(BfSizeArray const *sizeArray);
