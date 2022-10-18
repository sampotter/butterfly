#pragma once

#include "vec.h"

struct BfVecZero {
  BfVec super;
};

#define INTERFACE BF_INTERFACE_Vec
BF_DECLARE_INTERFACE(VecZero)
#undef INTERFACE

BfVec *bfVecZeroToVec(BfVecZero *vecZero);

BfVecZero const *bfVecConstToVecZeroConst(BfVec const *vec);

BfVecZero *bfVecZeroNew();
BfVecZero *bfVecZeroFromFile(char const *path, BfSize size);
void bfVecZeroInit(BfVecZero *vecZero, BfSize size);
void bfVecZeroDeinit(BfVecZero *vecZero);
void bfVecZeroDealloc(BfVecZero **vecZero);
void bfVecZeroDeinitAndDealloc(BfVecZero **vecZero);
