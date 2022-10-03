#pragma once

#include "mat.h"
#include "ptr_array.h"

typedef struct BfMatSum {
  BfMat super;
  BfPtrArray termArr;
} BfMatSum;

#define INTERFACE BF_INTERFACE_Mat
BF_DECLARE_INTERFACE(MatSum)
#undef INTERFACE

BfMat *bfMatSumToMat(BfMatSum *matSum);

BfMatSum const *bfMatConstToMatSumConst(BfMat const *mat);

BfMatSum *bfMatSumNew();
void bfMatSumInit(BfMatSum *sum);
void bfMatSumDeinit(BfMatSum *sum);
void bfMatSumDealloc(BfMatSum **sum);
void bfMatSumDeinitAndDealloc(BfMatSum **sum);
BfSize bfMatSumGetNumTerms(BfMatSum const *sum);
void bfMatSumAddTerm(BfMatSum *sum, BfMat *term);
BfMat *bfMatSumGetTerm(BfMatSum *sum, BfSize i);
BfMat const *bfMatSumGetTermConst(BfMatSum const *sum, BfSize i);
