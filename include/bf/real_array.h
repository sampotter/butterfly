#pragma once

#include "def.h"
#include "types.h"
#include "vec.h"

struct BfRealArray {
  BfReal *data;
  BfSize size;
  BfSize capacity;
};

BfRealArray *bfRealArrayNew();
BfRealArray *bfRealArrayNewWithDefaultCapacity();
void bfRealArrayInitWithDefaultCapacity(BfRealArray *realArray);
void bfRealArrayDeinit(BfRealArray *realArray);
void bfRealArrayDealloc(BfRealArray **realArray);
void bfRealArrayDeinitAndDealloc(BfRealArray **realArray);
void bfRealArrayExpandCapacity(BfRealArray *realArray, BfSize newCapacity);
void bfRealArrayAppend(BfRealArray *realArray, BfReal elt);
BfVec *bfRealArrayGetVecView(BfRealArray *realArray);
BfVec *bfRealArrayGetSubvecView(BfRealArray *realArray, BfSize i0, BfSize i1);
