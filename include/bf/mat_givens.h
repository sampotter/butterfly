#pragma once

#include "mat.h"

/** Interface: Mat */

void bfMatGivensComplexDelete(BfMatGivensComplex **matGivensComplex);
BfType bfMatGivensComplexGetType(BfMat const *mat);
BfSize bfMatGivensComplexGetNumRows(BfMat const *mat);
BfSize bfMatGivensComplexGetNumCols(BfMat const *mat);

struct BfMatGivensComplex {
  BfMat super;
  BfSize srcInd;
  BfSize elimInd;
  BfComplex c;
  BfComplex s;
};

BfMat *bfMatGivensComplexToMat(BfMatGivensComplex *mat);

BfMatGivensComplex const *bfMatConstToMatGivensComplexConst(BfMat const *mat);

BfMatGivensComplex *bfMatGivensComplexNew();
void bfMatGivensComplexInit(BfMatGivensComplex *mat, BfSize n, BfSize srcInd, BfSize elimInd, BfComplex c, BfComplex s);
void bfMatGivensComplexDeinit(BfMatGivensComplex *mat);
void bfMatGivensComplexDealloc(BfMatGivensComplex **mat);
void bfMatGivensComplexDeinitAndDealloc(BfMatGivensComplex **mat);
