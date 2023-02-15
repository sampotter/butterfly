#pragma once

#include "mat_block.h"

/** Interface: Mat */

BfMat *bfMatBlockDiagCopy(BfMat const *mat);
BfVec *bfMatBlockDiagGetRowCopy(BfMat const *mat, BfSize i);
void bfMatBlockDiagDelete(BfMat **mat);
BfType bfMatBlockDiagGetType(BfMat const *mat);
bool bfMatBlockDiagInstanceOf(BfMat const *mat, BfType type);
BfSize bfMatBlockDiagGetNumRows(BfMat const *mat);
BfSize bfMatBlockDiagGetNumCols(BfMat const *mat);
void bfMatBlockDiagScaleCols(BfMat *mat, BfVec const *vec);
BfMat *bfMatBlockDiagMul(BfMat const *mat, BfMat const *other);
void bfMatBlockDiagNegate(BfMat *mat);

/** Interface: MatBlock */

BfSize bfMatBlockDiagNumBlocks(BfMatBlockDiag const *matBlockDiag);
BfMat *bfMatBlockDiagGetBlockCopy(BfMatBlockDiag const *matBlockDiag, BfSize i, BfSize j);

struct BfMatBlockDiag {
  BfMatBlock super;
};

/* Upcasting: */
BfMat *bfMatBlockDiagToMat(BfMatBlockDiag *matBlockDiag);
BfMatBlock *bfMatBlockDiagToMatBlock(BfMatBlockDiag *);
BfMatBlock const *bfMatBlockDiagConstToMatBlockConst(BfMatBlockDiag const *);

/* Downcasting: */
BfMatBlockDiag *bfMatToMatBlockDiag(BfMat *mat);
BfMatBlockDiag const *bfMatConstToMatBlockDiagConst(BfMat const *mat);

BfMatBlockDiag *bfMatBlockDiagNew();
BfMatBlockDiag *bfMatBlockDiagNewFromBlocks(BfPtrArray *blocks);
void bfMatBlockDiagInit(BfMatBlockDiag *mat, BfSize numBlockRows, BfSize numBlockCols);
void bfMatBlockDiagDeinit(BfMatBlockDiag *mat);
void bfMatBlockDiagDealloc(BfMatBlockDiag **mat);
void bfMatBlockDiagDeinitAndDealloc(BfMatBlockDiag **mat);
BfMat *bfMatBlockDiagGetBlock(BfMatBlockDiag *matBlockDiag, BfSize i);
BfMat const *bfMatBlockDiagGetBlockConst(BfMatBlockDiag const *matBlockDiag, BfSize i);
void bfMatBlockDiagSetBlock(BfMatBlockDiag *matBlockDiag, BfSize i, BfMat *mat);
