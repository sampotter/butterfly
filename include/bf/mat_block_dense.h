#pragma once

#include "mat_block.h"

/** Interface: Mat */

BfMat *bfMatBlockDenseGetView(BfMatBlockDense *matBlockDense);
BfMat *bfMatBlockDenseCopy(BfMat const *mat);
BfMat *bfMatBlockDenseSteal(BfMatBlockDense *matBlockDense);
BfVec *bfMatBlockDenseGetRowCopy(BfMat const *mat, BfSize i);
void bfMatBlockDenseDelete(BfMat **mat);
BfType bfMatBlockDenseGetType(BfMat const *mat);
BfSize bfMatBlockDenseNumBytes(BfMat const *mat);
void bfMatBlockDenseSave(BfMatBlockDense const *matBlockDense, char const *path);
void bfMatBlockDenseDump(BfMatBlockDense const *matBlockDense, FILE *fp);
BfSize bfMatBlockDenseGetNumRows(BfMatBlockDense const *matBlockDense);
BfSize bfMatBlockDenseGetNumCols(BfMatBlockDense const *matBlockDense);
BfMat *bfMatBlockDenseGetRowRange(BfMatBlockDense *matBlockDense, BfSize i0, BfSize i1);
BfMat *bfMatBlockDenseGetRowRangeCopy(BfMatBlockDense const *matBlockDense, BfSize i0, BfSize i1);
void bfMatBlockDenseScaleCols(BfMatBlockDense *matBlockDense, BfVec const *vec);
void bfMatBlockDenseAddInplace(BfMatBlockDense *matBlockDense, BfMat const *otherMat);
BfMat *bfMatBlockDenseMul(BfMat const *mat, BfMat const *otherMat);
BfVec *bfMatBlockDenseMulVec(BfMatBlockDense const *matBlockDense, BfVec const *vec);
BfMat *bfMatBlockDenseRmul(BfMatBlockDense const *matBlockDense, BfMat const *mat);
BfVec *bfMatBlockDenseRmulVec(BfMatBlockDense const *matBlockDense, BfVec const *vec);
BfMat *bfMatBlockDenseToType(BfMatBlockDense const *matBlockDense, BfType type);
BfSizeArray *bfMatBlockDenseGetNonzeroColumnRanges(BfMatBlockDense const *matBlockDense);
void bfMatBlockDensePrintBlocksDeep(BfMatBlockDense const *matBlockDense, FILE *fp, BfSize i0, BfSize j0, BfSize depth);
void bfMatBlockDenseTranspose(BfMatBlockDense *matBlockDense);

/** Interface: MatBlock */

BfSize bfMatBlockDenseNumBlocks(BfMatBlockDense const *matBlockDense);
BfSize bfMatBlockDenseGetNumRowBlocks(BfMatBlockDense const *matBlockDense);
BfSize bfMatBlockDenseGetNumColBlocks(BfMatBlockDense const *matBlockDense);
BfSize bfMatBlockDenseGetRowOffset(BfMatBlockDense const *matBlockDense, BfSize i);
BfSize bfMatBlockDenseGetColOffset(BfMatBlockDense const *matBlockDense, BfSize j);
BfMat *bfMatBlockDenseGetBlock(BfMatBlockDense *matBlockDense, BfSize i, BfSize j);
BfMat const *bfMatBlockDenseGetBlockConst(BfMatBlockDense const *matBlockDense, BfSize i, BfSize j);
void bfMatBlockDenseGetRowSpan(BfMatBlockDense const *matBlockDense, BfSize i, BfSize *i0, BfSize *i1);
void bfMatBlockDenseGetColSpan(BfMatBlockDense const *matBlockDense, BfSize j, BfSize *j0, BfSize *j1);

/** Upcasting: MatBlockDense -> Mat */

BfMat *bfMatBlockDenseToMat(BfMatBlockDense *matBlockDense);
BfMat const *bfMatBlockDenseConstToMatConst(BfMatBlockDense const *matBlockDense);

/** Upcasting: MatBlockDense -> MatBlock */

BfMatBlock *bfMatBlockDenseToMatBlock(BfMatBlockDense *matBlockDense);
BfMatBlock const *bfMatBlockDenseConstToMatBlockConst(BfMatBlockDense const *matBlockDense);

/** Downcasting: Mat -> MatBlockDense */

BfMatBlockDense *bfMatToMatBlockDense(BfMat *mat);
BfMatBlockDense const *bfMatConstToMatBlockDenseConst(BfMat const *mat);

/** Implementation: MatBlockDense */

struct BfMatBlockDense {
  BfMatBlock super;
};

BfMatBlockDense *bfMatBlockDenseNew(void);
BfMatBlockDense *bfMatBlockDenseNewFromBlocks(BfSize numRowBlocks, BfSize numColBlocks, BfPtrArray *blocks, BfPolicy policy);
BfMatBlockDense *bfMatBlockDenseNewRowFromBlocks(BfPtrArray *blocks, BfPolicy policy);
BfMatBlockDense *bfMatBlockDenseNewColFromBlocks(BfPtrArray *blocks, BfPolicy policy);
void bfMatBlockDenseInit(BfMatBlockDense *mat, BfSize numBlockRows, BfSize numBlockCols);
void bfMatBlockDenseDeinit(BfMatBlockDense *mat);
void bfMatBlockDenseDealloc(BfMatBlockDense **mat);
void bfMatBlockDenseDeinitAndDealloc(BfMatBlockDense **mat);
void bfMatBlockDenseSetBlock(BfMatBlockDense *mat, BfSize i, BfSize j, BfMat *block);
void bfMatBlockDenseAppendColumn(BfMatBlockDense *mat, BfSize numBlockCols);
void bfMatBlockDenseAppendRow(BfMatBlockDense *mat, BfSize numBlockRows);
