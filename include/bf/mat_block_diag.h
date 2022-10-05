#pragma once

#include "mat_block.h"

struct BfMatBlockDiag {
  BfMatBlock super;
};

#define INTERFACE BF_INTERFACE_Mat
BF_DECLARE_INTERFACE(MatBlockDiag)
#undef INTERFACE

#define INTERFACE BF_INTERFACE_MatBlock
BF_DECLARE_INTERFACE(MatBlockDiag)
#undef INTERFACE

/* Upcasting: */
BfMat *bfMatBlockDiagToMat(BfMatBlockDiag *matBlockDiag);
BfMatBlock *bfMatBlockDiagToMatBlock(BfMatBlockDiag *);
BfMatBlock const *bfMatBlockDiagConstToMatBlockConst(BfMatBlockDiag const *);

/* Downcasting: */
BfMatBlockDiag *bfMatToMatBlockDiag(BfMat *mat);
BfMatBlockDiag const *bfMatConstToMatBlockDiagConst(BfMat const *mat);

BfMatBlockDiag *bfMatBlockDiagNew();
void bfMatBlockDiagInit(BfMatBlockDiag *mat, BfSize numBlockRows, BfSize numBlockCols);
void bfMatBlockDiagDeinit(BfMatBlockDiag *mat);
void bfMatBlockDiagDealloc(BfMatBlockDiag **mat);
void bfMatBlockDiagDeinitAndDealloc(BfMatBlockDiag **mat);
BfMat *bfMatBlockDiagGetBlock(BfMatBlockDiag *matBlockDiag, BfSize i);
BfMat const *bfMatBlockDiagGetBlockConst(BfMatBlockDiag const *matBlockDiag, BfSize i);
void bfMatBlockDiagSetBlock(BfMatBlockDiag *matBlockDiag, BfSize i, BfMat *mat);
