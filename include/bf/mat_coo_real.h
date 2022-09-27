#pragma once

#include "mat.h"

struct BfMatCooReal {
  BfMat super;
  BfSize numElts;
  BfSize *rowInd;
  BfSize *colInd;
  BfReal *value;
};

#define INTERFACE BF_INTERFACE_Mat
BF_DECLARE_INTERFACE(MatCooReal)
#undef INTERFACE

/* Upcasting: */
BfMat *bfMatCooRealToMat(BfMatCooReal *matCooReal);
BfMat const *bfMatCooRealConstToMatConst(BfMatCooReal const *matCooReal);

/* Downcasting: */
BfMatCooReal *bfMatToMatCooReal(BfMat *mat);
BfMatCooReal const *bfMatConstToMatCooRealConst(BfMat const *mat);

BfMatCooReal *bfMatCooRealNew();
void bfMatCooRealInitEmpty(BfMatCooReal *mat, BfSize numRows, BfSize numCols, BfSize numElts);
void bfMatCooRealDeinit(BfMatCooReal *mat);
void bfMatCooRealDealloc(BfMatCooReal **mat);
void bfMatCooRealDeinitAndDealloc(BfMatCooReal **mat);
