#pragma once

#include "mat.h"

/** Interface: MatCooReal */

void bfMatCooRealDelete(BfMat **mat);
BfType bfMatCooRealGetType(BfMat const *mat);
BfSize bfMatCooRealGetNumRows(BfMat const *mat);
BfSize bfMatCooRealGetNumCols(BfMat const *mat);
BfMat *bfMatCooRealGetRowRangeCopy(BfMat const *mat, BfSize i0, BfSize i1);
BfMat *bfMatCooRealGetColRangeCopy(BfMat const *mat, BfSize j0, BfSize j1);
void bfMatCooRealPermuteRows(BfMat *mat, BfPerm const *perm);
void bfMatCooRealPermuteCols(BfMat *mat, BfPerm const *perm);
bool bfMatCooRealIsZero(BfMat const *mat);

struct BfMatCooReal {
  BfMat super;
  BfSize numElts;
  BfSize *rowInd;
  BfSize *colInd;
  BfReal *value;
};

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
