#pragma once

#include "vec.h"

void bfVecZeroDelete(BfVec **mat);
BfType bfVecZeroGetType(BfVec const *vec);
BfVec *bfVecZeroConcat(BfVec const *vec, BfVec const *otherVec);

struct BfVecZero {
  BfVec super;
};

BfVec *bfVecZeroToVec(BfVecZero *vecZero);

BfVecZero const *bfVecConstToVecZeroConst(BfVec const *vec);

BfVecZero *bfVecZeroNew(void);
BfVecZero *bfVecZeroFromFile(char const *path, BfSize size);
void bfVecZeroInit(BfVecZero *vecZero, BfSize size);
void bfVecZeroDeinit(BfVecZero *vecZero);
void bfVecZeroDealloc(BfVecZero **vecZero);
void bfVecZeroDeinitAndDealloc(BfVecZero **vecZero);
