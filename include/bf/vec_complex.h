#pragma once

#include "vec.h"

struct BfVecComplex {
  BfVec super;
  BfSize stride;
  BfComplex *data;
};

#define INTERFACE BF_INTERFACE_Vec
BF_DECLARE_INTERFACE(VecComplex)
#undef INTERFACE

BfVec *bfVecComplexToVec(BfVecComplex *vecComplex);

BfVecComplex *bfVecToVecComplex(BfVec *vec);
BfVecComplex const *bfVecConstToVecComplexConst(BfVec const *vec);

BfVecComplex *bfVecComplexNew();
BfVecComplex *bfVecComplexFromFile(char const *path, BfSize size);
void bfVecComplexInit(BfVecComplex *vecComplex, BfSize size);
void bfVecComplexInitView(BfVecComplex *vecComplex, BfSize size, BfSize stride, BfComplex *data);
void bfVecComplexDeinit(BfVecComplex *vecComplex);
void bfVecComplexDealloc(BfVecComplex **vecComplex);
void bfVecComplexDeinitAndDealloc(BfVecComplex **vecComplex);
