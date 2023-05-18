#pragma once

#include "mat.h"

/** Interface: MatZero */

BfVec *bfMatZeroGetRowCopy(BfMat const *mat, BfSize i);
void bfMatZeroDelete(BfMatZero **matZero);
BfType bfMatZeroGetType(BfMat const *mat);
BfSize bfMatZeroGetNumRows(BfMat const *mat);
BfSize bfMatZeroGetNumCols(BfMat const *mat);

struct BfMatZero {
  BfMat super;
};

BfMat *bfMatZeroToMat(BfMatZero *mat);
BfMat const *bfMatZeroConstToMatConst(BfMatZero const *mat);

BfMatZero const *bfMatConstToMatZeroConst(BfMat const *mat);

BfMatZero *bfMatZeroNew();
void bfMatZeroInit(BfMatZero *mat, BfSize numRows, BfSize numCols);
void bfMatZeroDeinit(BfMatZero *mat);
void bfMatZeroDealloc(BfMatZero **mat);
void bfMatZeroDeinitAndDealloc(BfMatZero **mat);
