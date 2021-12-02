#pragma once

#include "def.h"
#include "error.h"

static BfSize const BF_ARRAY_DEFAULT_CAPACITY = 1024;

typedef struct BfPtrArray {
  BfPtr *data;
  BfSize capacity, num_elts;
} BfPtrArray;

typedef enum BfError (*BfPtrFunc)(BfPtr elt_ptr, BfPtr arg_ptr);

enum BfError
bfInitPtrArray(BfPtrArray *arr, BfSize capacity);

enum BfError
bfInitPtrArrayWithDefaultCapacity(BfPtrArray *arr);

enum BfError
bfFreePtrArray(BfPtrArray *arr);

BfSize bfPtrArraySize(BfPtrArray const *arr);

enum BfError
bfPtrArrayAppend(BfPtrArray *arr, BfPtr ptr);

enum BfError
bfPtrArrayGet(BfPtrArray const *arr, BfSize pos, BfPtr *ptr);

enum BfError
bfPtrArrayGetFirst(BfPtrArray const *arr, BfPtr *ptr);

enum BfError
bfPtrArrayGetLast(BfPtrArray const *arr, BfPtr *ptr);

enum BfError
bfMapPtrArray(BfPtrArray *arr, BfPtrFunc func, BfPtr ptr);
