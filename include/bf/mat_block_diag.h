#pragma once

#include "mat_block.h"

/** Interface: Mat */

BfMat *bfMatBlockDiagGetView(BfMatBlockDiag *matBlockDiag);
BfMat *bfMatBlockDiagCopy(BfMat const *mat);
BfMat *bfMatBlockDiagSteal(BfMatBlockDiag *matBlockDiag);
BfVec *bfMatBlockDiagGetRowCopy(BfMat const *mat, BfSize i);
void bfMatBlockDiagDelete(BfMat **mat);
BfType bfMatBlockDiagGetType(BfMat const *mat);
BfSize bfMatBlockDiagNumBytes(BfMatBlockDiag const *matBlockDiag);
void bfMatBlockDiagDump(BfMatBlockDiag const *matBlockDiag, FILE *fp);
bool bfMatBlockDiagInstanceOf(BfMat const *mat, BfType type);
BfSize bfMatBlockDiagGetNumRows(BfMatBlockDiag const *matBlockDiag);
BfSize bfMatBlockDiagGetNumCols(BfMatBlockDiag const *matBlockDiag);
BfMat *bfMatBlockDiagGetRowRangeCopy(BfMatBlockDiag const *matBlockDiag, BfSize i0, BfSize i1);
void bfMatBlockDiagScaleCols(BfMat *mat, BfVec const *vec);
BfMat *bfMatBlockDiagMul(BfMat const *mat, BfMat const *other);
BfVec *bfMatBlockDiagMulVec(BfMatBlockDiag const *matBlockDiag, BfVec const *vec);
BfVec *bfMatBlockDiagRmulVec(BfMatBlockDiag const *matBlockDiag, BfVec const *vec);
void bfMatBlockDiagNegate(BfMat *mat);
void bfMatBlockDiagPrintBlocksDeep(BfMatBlockDiag const *matBlockDiag, FILE *fp, BfSize i0, BfSize j0, BfSize depth);
BfMat *bfMatBlockDiagSolve(BfMatBlockDiag const *matBlockDiag, BfMat const *mat);
void bfMatBlockDiagTranspose(BfMatBlockDiag *matBlockDiag);

/** Interface: MatBlock */

BfSize bfMatBlockDiagNumBlocks(BfMatBlockDiag const *matBlockDiag);
BfSize bfMatBlockDiagGetNumRowBlocks(BfMatBlockDiag const *matBlockDiag);
BfSize bfMatBlockDiagGetNumColBlocks(BfMatBlockDiag const *matBlockDiag);
BfMat *bfMatBlockDiagGetBlockCopy(BfMatBlockDiag const *matBlockDiag, BfSize i, BfSize j);

struct BfMatBlockDiag {
  BfMatBlock super;
};

/** Upcasting: MatBlockDiag -> Mat */

BfMat *bfMatBlockDiagToMat(BfMatBlockDiag *matBlockDiag);
BfMat const *bfMatBlockDiagConstToMatConst(BfMatBlockDiag const *matBlockDiag);

/** Upcasting: MatBlockDiag -> MatBlock */

BfMatBlock *bfMatBlockDiagToMatBlock(BfMatBlockDiag *);
BfMatBlock const *bfMatBlockDiagConstToMatBlockConst(BfMatBlockDiag const *);

/** Downcasting: Mat -> MatBlockDiag */

BfMatBlockDiag *bfMatToMatBlockDiag(BfMat *mat);
BfMatBlockDiag const *bfMatConstToMatBlockDiagConst(BfMat const *mat);

/** Implementation: MatBlockDiag */

BfMatBlockDiag *bfMatBlockDiagNew(void);
BfMatBlockDiag *bfMatBlockDiagNewFromBlocks(BfPtrArray *blocks, BfPolicy policy);
void bfMatBlockDiagInit(BfMatBlockDiag *mat, BfSize numBlockRows, BfSize numBlockCols);
void bfMatBlockDiagDeinit(BfMatBlockDiag *mat);
void bfMatBlockDiagDealloc(BfMatBlockDiag **mat);
void bfMatBlockDiagDeinitAndDealloc(BfMatBlockDiag **mat);
BfMat *bfMatBlockDiagGetBlock(BfMatBlockDiag *matBlockDiag, BfSize i);
BfMat const *bfMatBlockDiagGetBlockConst(BfMatBlockDiag const *matBlockDiag, BfSize i);
void bfMatBlockDiagSetBlock(BfMatBlockDiag *matBlockDiag, BfSize i, BfMat *mat);
