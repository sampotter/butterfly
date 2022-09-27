#pragma once

#include "vec.h"

struct BfVecReal {
  BfVec super;
  BfSize stride;
  BfReal *data;
};

#define INTERFACE BF_INTERFACE_Vec
BF_DECLARE_INTERFACE(VecReal)
#undef INTERFACE

BfVec *bfVecRealToVec(BfVecReal *vecReal);

BfVecReal *bfVecToVecReal(BfVec *vec);
BfVecReal const *bfVecConstToVecRealConst(BfVec const *vec);

BfVecReal *bfVecRealNew();
BfVecReal *bfVecRealFromFile(char const *path, BfSize size);
void bfVecRealInit(BfVecReal *vecReal, BfSize size);
void bfVecRealInitView(BfVecReal *vecReal, BfSize size, BfSize stride, BfReal *data);
void bfVecRealDeinit(BfVecReal *vecReal);
void bfVecRealDealloc(BfVecReal **vecReal);
void bfVecRealDeinitAndDealloc(BfVecReal **vecReal);
