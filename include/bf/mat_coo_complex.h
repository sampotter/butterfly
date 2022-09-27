#pragma once

#include "mat.h"

struct BfMatCooComplex {
  BfMat super;
  BfSize numElts;
  BfSize *rowInd;
  BfSize *colInd;
  BfComplex *value;
};

#define INTERFACE BF_INTERFACE_Mat
BF_DECLARE_INTERFACE(MatCooComplex)
#undef INTERFACE

/* Upcasting: */
BfMat *bfMatCooComplexToMat(BfMatCooComplex *matCooComplex);
BfMat const *bfMatCooComplexConstToMatConst(BfMatCooComplex const *matCooComplex);

/* Downcasting: */
BfMatCooComplex *bfMatToMatCooComplex(BfMat *mat);
BfMatCooComplex const *bfMatConstToMatCooComplexConst(BfMat const *mat);

BfMatCooComplex *bfMatCooComplexNew();
void bfMatCooComplexInitEmpty(BfMatCooComplex *mat, BfSize numRows, BfSize numCols, BfSize numElts);
void bfMatCooComplexDeinit(BfMatCooComplex *mat);
void bfMatCooComplexDealloc(BfMatCooComplex **mat);
void bfMatCooComplexDeinitAndDealloc(BfMatCooComplex **mat);
