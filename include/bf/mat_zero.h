#pragma once

#include "mat.h"

struct BfMatZero {
  BfMat super;
};

#define INTERFACE BF_INTERFACE_Mat
BF_DECLARE_INTERFACE(MatZero)
#undef INTERFACE

BfMatZero *bfMatZeroNew();
void bfMatZeroInit(BfMatZero *mat, BfSize numRows, BfSize numCols);
void bfMatZeroDeinit(BfMatZero *mat);
void bfMatZeroDealloc(BfMatZero **mat);
void bfMatZeroDeinitAndDealloc(BfMatZero **mat);
