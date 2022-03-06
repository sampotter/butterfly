#pragma once

#include "mat_block.h"

struct BfMatBlockDiag {
  BfMatBlock super;
};

BfMatBlockDiag *bfMatBlockDiagNew();
void bfMatBlockDiagInit(BfMatBlockDiag *mat, BfSize numBlockRows, BfSize numBlockCols);
BfMat *bfMatBlockDiagGetMatPtr(BfMatBlockDiag *mat);

/* BfMat interface: */
void bfMatBlockDiagDeinit(BfMatBlockDiag *mat);
void bfMatBlockDiagDelete(BfMatBlockDiag **mat);
void bfMatBlockDiagDeinitAndDelete(BfMatBlockDiag **mat);
BfMatType bfMatBlockDiagGetType(BfMatBlockDiag const *mat);
BfSize bfMatBlockDiagNumBytes(BfMatBlockDiag const *mat);
void bfMatBlockDiagSave(BfMatBlockDiag const *mat, char const *path);
BfMat *bfMatBlockDiagMul(BfMatBlockDiag const *op1, BfMat const *op2);
BfMat *bfMatBlockDiagLstSq(BfMatBlockDiag const *op1, BfMat const *op2);

/* BfMatBlock interface: */
BfSize bfMatBlockDiagNumBlocks(BfMatBlockDiag const *mat);
