#pragma once

#include "mat.h"

struct BfMatGivensComplex {
  BfMat super;
  BfSize srcInd;
  BfSize elimInd;
  BfReal c;
  BfComplex s;
};

#define INTERFACE BF_INTERFACE_Mat
BF_DECLARE_INTERFACE(MatGivensComplex)
#undef INTERFACE

BfMat *bfMatGivensComplexToMat(BfMatGivensComplex *mat);

BfMatGivensComplex const *bfMatConstToMatGivensComplexConst(BfMat const *mat);

BfMatGivensComplex *bfMatGivensComplexNew();
void bfMatGivensComplexInit(BfMatGivensComplex *mat, BfSize n, BfSize srcInd, BfSize elimInd, BfReal c, BfComplex s);
void bfMatGivensComplexDeinit(BfMatGivensComplex *mat);
void bfMatGivensComplexDealloc(BfMatGivensComplex **mat);
void bfMatGivensComplexDeinitAndDealloc(BfMatGivensComplex **mat);
