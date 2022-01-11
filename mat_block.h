#pragma once

#include "mat.h"

typedef struct BfMatBlock {
  BfMat super;

  BfSize numBlockRows;
  BfSize numBlockCols;
} BfMatBlock;

BfSize bfMatBlockGetNumBlockRows(BfMatBlock const *mat);
