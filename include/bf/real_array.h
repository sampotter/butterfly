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

BfRealArray *bfRealArrayNew(void);
BfRealArray *bfRealArrayNewFromVecReal(BfVecReal const *vecReal, BfPolicy policy);
BfRealArray *bfRealArrayNewWithCapacity(BfSize capacity);
BfRealArray *bfRealArrayNewWithDefaultCapacity(void);
BfRealArray *bfRealArrayNewWithValue(BfSize size, BfReal value);
BfRealArray *bfRealArrayNewFromFile(char const *path);
BfRealArray *bfRealArrayCopy(BfRealArray const *realArray);
void bfRealArrayInitWithCapacity(BfRealArray *realArray, BfSize capacity);
void bfRealArrayInitWithDefaultCapacity(BfRealArray *realArray);
void bfRealArrayInitWithValue(BfRealArray *realArray, BfSize size, BfReal value);
void bfRealArrayInitCopy(BfRealArray *realArray, BfRealArray const *otherRealArray);
void bfRealArrayDeinit(BfRealArray *realArray);
void bfRealArrayDealloc(BfRealArray **realArray);
void bfRealArrayDeinitAndDealloc(BfRealArray **realArray);
void bfRealArrayExpandCapacity(BfRealArray *realArray, BfSize newCapacity);
void bfRealArrayShrinkCapacityToSize(BfRealArray *realArray);
void bfRealArrayAppend(BfRealArray *realArray, BfReal elt);
void bfRealArrayExtend(BfRealArray *realArray, BfRealArray const *otherRealArray);
BfVec *bfRealArrayGetVecView(BfRealArray *realArray);
BfVec *bfRealArrayGetSubvecView(BfRealArray *realArray, BfSize i0, BfSize i1);
BfReal bfRealArrayGetValue(BfRealArray const *realArray, BfSize i);
void bfRealArraySetValue(BfRealArray *realArray, BfSize i, BfReal value);
void bfRealArrayGetValues(BfRealArray const *realArray, BfSize n, BfSize const *inds, BfReal *values);
void bfRealArrayInsert(BfRealArray *realArray, BfSize i, BfReal value);
BfSize bfRealArrayGetSize(BfRealArray const *realArray);
bool bfRealArrayIsEmpty(BfRealArray const *realArray);
void bfRealArraySave(BfRealArray const *realArray, char const *path);
void bfRealArrayNegate(BfRealArray *realArray);
void bfRealArrayPermute(BfRealArray *realArray, BfPerm const *perm);
BfReal *bfRealArrayGetDataPtr(BfRealArray *realArray);
BfReal *bfRealArrayStealPtr(BfRealArray *realArray);
BfPerm *bfRealArrayArgsort(BfRealArray const *realArray);
BfRealArray *bfRealArrayGetRangeView(BfRealArray const *realArray, BfSize i0, BfSize i1);
