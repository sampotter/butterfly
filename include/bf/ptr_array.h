#pragma once

#include "def.h"
#include "error.h"
#include "types.h"

struct BfPtrArray {
  BfPtr *data;
  BfSize capacity;
  BfSize num_elts;
  bool isView;
};

typedef void (*BfPtrFunc)(BfPtr elt_ptr, BfPtr arg_ptr);
typedef int (*BfPtrCmp)(BfPtr, BfPtr);

BfPtrArray *bfPtrArrayNewWithDefaultCapacity(void);
BfPtrArray bfGetUninitializedPtrArray(void);
void bfInitPtrArray(BfPtrArray *arr, BfSize capacity);
void bfInitPtrArrayWithDefaultCapacity(BfPtrArray *arr);
void bfMakeEmptyPtrArrayView(BfPtrArray *arr);
void bfPtrArrayDeinit(BfPtrArray *arr);
void bfPtrArrayDealloc(BfPtrArray **arr);
void bfPtrArrayDelete(BfPtrArray **arr);
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
void bfPtrArraySort(BfPtrArray *arr, BfPtrCmp ptrCmp, BfPtr aux);
void bfPtrArrayReverse(BfPtrArray *arr);
void bfPtrArrayCopyData(BfPtrArray const *ptrArray, BfPtr *dst);
void bfPtrArrayRemove(BfPtrArray *ptrArray, BfSize i);

struct BfConstPtrArray {
  BfConstPtr *data;
  BfSize capacity, num_elts;
  bool isView;
};

BfConstPtrArray *bfConstPtrArrayNewWithDefaultCapacity(void);
void bfConstPtrArrayInit(BfConstPtrArray *arr, BfSize capacity);
void bfConstPtrArrayInitWithDefaultCapacity(BfConstPtrArray *arr);
void bfConstPtrArrayDeinit(BfConstPtrArray *arr);
void bfConstPtrArrayDealloc(BfConstPtrArray **arr);
void bfConstPtrArrayDeinitAndDealloc(BfConstPtrArray **arr);
bool bfConstPtrArrayIsEmpty(BfConstPtrArray const *arr);
BfSize bfConstPtrArraySize(BfConstPtrArray const *arr);
void bfConstPtrArrayAppend(BfConstPtrArray *arr, BfConstPtr ptr);
void bfConstPtrArrayExtend(BfConstPtrArray *arr, BfConstPtrArray const *otherArr);
BfConstPtr bfConstPtrArrayGet(BfConstPtrArray const *arr, BfSize pos);
BfConstPtr bfConstPtrArrayPopLast(BfConstPtrArray *arr);
BfConstPtr bfConstPtrArrayGetFirst(BfConstPtrArray const *arr);
BfConstPtr bfConstPtrArrayGetLast(BfConstPtrArray const *arr);
