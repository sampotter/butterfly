#pragma once

#include "fac.h"

struct BfFacSpan {
  BfSize numFacs;
  BfFac **fac;
};

BfFacSpan *bfFacSpanNewFromPtrArray(BfPtrArray const *ptrArray);
void bfFacSpanInitFromPtrArray(BfFacSpan *facSpan, BfPtrArray const *ptrArray);
void bfFacSpanDeinit(BfFacSpan *facSpan);
void bfFacSpanDealloc(BfFacSpan **facSpan);
void bfFacSpanDelete(BfFacSpan **facSpan);
BfMat *bfFacSpanGetMat(BfFacSpan const *facSpan, BfPolicy policy);
