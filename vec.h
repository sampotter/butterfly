#pragma once

#include "def.h"
#include "dtype.h"

enum BfVecProps {
  BF_VEC_PROP_NONE = 0,
  BF_VEC_PROP_VIEW = (1 << 0)
};

typedef struct BfVec {
  enum BfDtypes dtype;
  enum BfVecProps props;
  BfSize size;
  BfSize stride;
  BfPtr data;
} BfVec;

BfPtr bfVecGetEltPtr(BfVec *x, BfSize i);
BfReal bfVecDist(BfVec const *x, BfVec const *y);
void bfVecScaleByReal(BfVec *x, BfReal scale);
