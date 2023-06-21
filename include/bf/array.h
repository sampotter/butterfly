#pragma once

#include "def.h"
#include "types.h"

BfArray *bfArrayCopy(BfArray const *array);
BfArray *bfArrayNewUninitialized(void);
BfArray *bfArrayNewEmpty(BfSize eltSize);
BfArray *bfArrayNewWithValue(BfSize eltSize, BfSize size, BfConstPtr eltPtr);
void bfArrayInitEmpty(BfArray *array, BfSize eltSize);
void bfArrayInitWithValue(BfArray *array, BfSize eltSize, BfSize size, BfConstPtr eltPtr);
void bfArrayDeinit(BfArray *array);
void bfArrayDealloc(BfArray **array);
void bfArrayDeinitAndDealloc(BfArray **array);
BfSize bfArrayGetSize(BfArray const *array);
void bfArrayGet(BfArray const *array, BfSize i, BfPtr eltPtr);
BfPtr bfArrayGetPtr(BfArray *array, BfSize i);
BfSize bfArrayFindSorted(BfArray const *array, BfConstPtr eltPtr, BfCompar compar);
void bfArrayInsert(BfArray *array, BfSize i, BfConstPtr eltPtr);
void bfArraySave(BfArray const *array, char const *path);
bool bfArrayIsEmpty(BfArray const *array);
void bfArrayRemove(BfArray *array, BfSize i);
