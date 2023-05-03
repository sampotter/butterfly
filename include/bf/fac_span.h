#pragma once

#include "fac.h"

struct BfFacSpan {
  BfSize numFacs;
  BfFac **fac;
};

BfFacSpan *bfFacSpanNewFromPtrArray(BfPtrArray const *ptrArray);
void bfFacSpanInitFromPtrArray(BfFacSpan *facSpan, BfPtrArray const *ptrArray);
BfMat *bfFacSpanGetMat(BfFacSpan const *facSpan);
