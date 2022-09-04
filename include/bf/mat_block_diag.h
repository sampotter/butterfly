#pragma once

#include "mat_block.h"

struct BfMatBlockDiag {
  BfMatBlock super;
};

#define INTERFACE BF_INTERFACE_Mat
BF_DECLARE_INTERFACE(MatBlockDiag);
#undef INTERFACE

#define INTERFACE BF_INTERFACE_MatBlock
BF_DECLARE_INTERFACE(MatBlockDiag);
#undef INTERFACE

/* Upcasting: */
BfMat *bfMatBlockDiagToMat(BfMatBlockDiag *matBlockDiag);

/* Downcasting: */
BfMatBlockDiag const *bfMatConstToMatBlockDiagConst(BfMat const *mat);

BfMatBlockDiag *bfMatBlockDiagNew();
void bfMatBlockDiagInit(BfMatBlockDiag *mat, BfSize numBlockRows, BfSize numBlockCols);
void bfMatBlockDiagDeinit(BfMatBlockDiag *mat);
void bfMatBlockDiagDealloc(BfMatBlockDiag **mat);
void bfMatBlockDiagDeinitAndDealloc(BfMatBlockDiag **mat);
