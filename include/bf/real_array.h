#pragma once

#include "def.h"
#include "types.h"
#include "vec.h"

struct BfRealArray {
  BfReal *data;
  BfSize size;
  BfSize capacity;
  bool isView;
};

BfRealArray *bfRealArrayNew();
BfRealArray *bfRealArrayNewFromVecReal(BfVecReal const *vecReal, BfPolicy policy);
BfRealArray *bfRealArrayNewWithDefaultCapacity();
BfRealArray *bfRealArrayNewWithValue(BfSize size, BfReal value);
void bfRealArrayInitWithDefaultCapacity(BfRealArray *realArray);
void bfRealArrayInitWithValue(BfRealArray *realArray, BfSize size, BfReal value);
void bfRealArrayDeinit(BfRealArray *realArray);
void bfRealArrayDealloc(BfRealArray **realArray);
void bfRealArrayDeinitAndDealloc(BfRealArray **realArray);
void bfRealArrayExpandCapacity(BfRealArray *realArray, BfSize newCapacity);
void bfRealArrayAppend(BfRealArray *realArray, BfReal elt);
BfVec *bfRealArrayGetVecView(BfRealArray *realArray);
BfVec *bfRealArrayGetSubvecView(BfRealArray *realArray, BfSize i0, BfSize i1);
BfReal bfRealArrayGetValue(BfRealArray const *realArray, BfSize i);
void bfRealArrayGetValues(BfRealArray const *realArray, BfSize n, BfSize const *inds, BfReal *values);
void bfRealArrayInsert(BfRealArray *realArray, BfSize i, BfReal value);
BfSize bfRealArrayGetSize(BfRealArray const *realArray);
