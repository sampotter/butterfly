#pragma once

#include "mat.h"

/** Interface: MatDiff */

BfType bfMatDiffGetType(BfMatDiff const *matDiff);
BfSize bfMatDiffGetNumRows(BfMatDiff const *matDiff);
BfSize bfMatDiffGetNumCols(BfMatDiff const *matDiff);
BfMat *bfMatDiffMul(BfMatDiff const *matDiff, BfMat const *otherMat);
BfMat *bfMatDiffToType(BfMatDiff const *matDiff, BfType type);

/** Upcasting: MatDiff -> Mat */

BfMat *bfMatDiffToMat(BfMatDiff *matDiff);
BfMat const *bfMatDiffConstToMatConst(BfMatDiff const  *matDiff);

/** Downcasting: Mat -> MatDiff */

BfMatDiff const *bfMatConstToMatDiffConst(BfMat const *mat);

/** Implementation: MatDiff */

struct BfMatDiff {
  BfMat super;
  BfMat *first;
  BfMat *second;
};

BfMatDiff *bfMatDiffAlloc(void);
BfMatDiff *bfMatDiffNew(BfMat *first, BfMat *second, BfPolicy policy);
void bfMatDiffInit(BfMatDiff *diff, BfMat *first, BfMat *second, BfPolicy policy);
void bfMatDiffDeinit(BfMatDiff *diff);
void bfMatDiffDealloc(BfMatDiff **diff);
void bfMatDiffDeinitAndDealloc(BfMatDiff **diff);
BfMat *bfMatDiffGetFirst(BfMatDiff *matDiff);
BfMat *bfMatDiffGetSecond(BfMatDiff *matDiff);
