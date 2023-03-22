#pragma once

#include "mat.h"

/** Interface: MatCooComplex */

void bfMatCooComplexDelete(BfMat **mat);
BfType bfMatCooComplexGetType(BfMat const *mat);
BfSize bfMatCooComplexGetNumRows(BfMat const *mat);
BfSize bfMatCooComplexGetNumCols(BfMat const *mat);
BfMat *bfMatCooComplexGetRowRangeCopy(BfMat const *mat, BfSize i0, BfSize i1);
BfMat *bfMatCooComplexGetColRangeCopy(BfMat const *mat, BfSize j0, BfSize j1);
void bfMatCooComplexPermuteRows(BfMat *mat, BfPerm const *perm);
void bfMatCooComplexPermuteCols(BfMat *mat, BfPerm const *perm);
BfMat *bfMatCooComplexMul(BfMat const *mat, BfMat const *otherMat);
bool bfMatCooComplexIsZero(BfMat const *mat);

struct BfMatCooComplex {
  BfMat super;
  BfSize numElts;
  BfSize *rowInd;
  BfSize *colInd;
  BfComplex *value;
};

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
