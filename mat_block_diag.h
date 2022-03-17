#pragma once

#include "mat_block.h"

struct BfMatBlockDiag {
  BfMatBlock super;
};

BF_DECLARE_INTERFACE_MAT(MatBlockDiag);

BfMatBlockDiag *bfMatBlockDiagNew();
void bfMatBlockDiagInit(BfMatBlockDiag *mat, BfSize numBlockRows, BfSize numBlockCols);
BfMat *bfMatBlockDiagGetMatPtr(BfMatBlockDiag *mat);
BfMat const *bfMatBlockDiagGetMatConstPtr(BfMatBlockDiag const *mat);

/* BfMatBlock interface: */
BfSize bfMatBlockDiagNumBlocks(BfMatBlockDiag const *mat);
BfMat const *bfMatBlockDiagGetBlock(BfMatBlockDiag const *mat, BfSize i, BfSize j);
