#pragma once

#include "def.h"
#include "error.h"
#include "types.h"

static BfSize const BF_ARRAY_DEFAULT_CAPACITY = 1024;

enum BfPtrArrayFlags {
  BF_PTR_ARRAY_FLAG_NONE = 0,
  BF_PTR_ARRAY_FLAG_VIEW = 1 << 0
};

struct BfPtrArray {
  enum BfPtrArrayFlags flags;
  BfPtr *data;
  BfSize capacity, num_elts;
};

typedef void (*BfPtrFunc)(BfPtr elt_ptr, BfPtr arg_ptr);
typedef int (*BfPtrCmp)(BfPtr, BfPtr);

BfPtrArray bfGetUninitializedPtrArray();
void bfInitPtrArray(BfPtrArray *arr, BfSize capacity);
void bfInitPtrArrayWithDefaultCapacity(BfPtrArray *arr);
void bfMakeEmptyPtrArrayView(BfPtrArray *arr);
void bfPtrArrayDeinit(BfPtrArray *arr);
BfPtrArray bfPtrArrayCopy(BfPtrArray *arr);
BfSize bfPtrArraySize(BfPtrArray const *arr);
bool bfPtrArrayIsEmpty(BfPtrArray const *arr);
void bfPtrArrayAppend(BfPtrArray *arr, BfPtr ptr);
BfPtr bfPtrArrayGet(BfPtrArray const *arr, BfSize pos);
void bfPtrArrayGetFirst(BfPtrArray const *arr, BfPtr *ptr);
void bfPtrArrayGetLast(BfPtrArray const *arr, BfPtr *ptr);
void bfPtrArrayGetRangeView(BfPtrArray const *arr, BfSize start, BfSize end, BfPtrArray *view);
void bfMapPtrArray(BfPtrArray *arr, BfPtrFunc func, BfPtr ptr);
BfPtr bfPtrArrayPopLast(BfPtrArray *arr);
void bfPtrArraySort(BfPtrArray *arr, BfPtrCmp ptrCmp);
void bfPtrArrayReverse(BfPtrArray *arr);
void bfPtrArrayCopyData(BfPtrArray const *ptrArray, BfPtr *dst);

struct BfConstPtrArray {
  enum BfPtrArrayFlags flags;
  BfConstPtr *data;
  BfSize capacity, num_elts;
};

void bfConstPtrArrayDeinit(BfConstPtrArray *arr);
bool bfConstPtrArrayIsEmpty(BfConstPtrArray const *arr);
void bfConstPtrArrayAppend(BfConstPtrArray *arr, BfConstPtr ptr);
BfConstPtr bfConstPtrArrayPopLast(BfConstPtrArray *arr);
