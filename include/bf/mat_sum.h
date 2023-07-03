#pragma once

#include "mat.h"
#include "ptr_array.h"

/** Interface: MatSum */

BfType bfMatSumGetType(BfMat const *mat);
BfSize bfMatSumGetNumRows(BfMat const *mat);
BfSize bfMatSumGetNumCols(BfMat const *mat);
BfMat *bfMatSumMul(BfMat const *mat, BfMat const *otherMat);

/** Implementation: MatSum */

typedef struct BfMatSum {
  BfMat super;
  BfPtrArray termArr;
} BfMatSum;

BfMat *bfMatSumToMat(BfMatSum *matSum);

BfMatSum const *bfMatConstToMatSumConst(BfMat const *mat);

BfMatSum *bfMatSumNew(void);
void bfMatSumInit(BfMatSum *sum);
void bfMatSumDeinit(BfMatSum *sum);
void bfMatSumDealloc(BfMatSum **sum);
void bfMatSumDeinitAndDealloc(BfMatSum **sum);
BfSize bfMatSumGetNumTerms(BfMatSum const *sum);
void bfMatSumAddTerm(BfMatSum *sum, BfMat *term);
BfMat *bfMatSumGetTerm(BfMatSum *sum, BfSize i);
BfMat const *bfMatSumGetTermConst(BfMatSum const *sum, BfSize i);
