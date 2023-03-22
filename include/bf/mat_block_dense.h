#pragma once

#include "mat_block.h"

/** Interface: Mat */

BfMat *bfMatBlockDenseCopy(BfMat const *mat);
BfVec *bfMatBlockDenseGetRowCopy(BfMat const *mat, BfSize i);
void bfMatBlockDenseDelete(BfMat **mat);
BfType bfMatBlockDenseGetType(BfMat const *mat);
BfSize bfMatBlockDenseNumBytes(BfMat const *mat);
BfSize bfMatBlockDenseGetNumRows(BfMat const *mat);
BfSize bfMatBlockDenseGetNumCols(BfMat const *mat);
BfMat *bfMatBlockDenseGetRowRange(BfMat *mat, BfSize i0, BfSize i1);
BfMat *bfMatBlockDenseGetRowRangeCopy(BfMatBlockDense const *matBlockDense, BfSize i0, BfSize i1);
void bfMatBlockDenseScaleCols(BfMatBlockDense *matBlockDense, BfVec const *vec);
void bfMatBlockDenseAddInplace(BfMatBlockDense *matBlockDense, BfMat const *otherMat);
BfMat *bfMatBlockDenseMul(BfMat const *mat, BfMat const *otherMat);
BfVec *bfMatBlockDenseMulVec(BfMatBlockDense const *matBlockDense, BfVec const *vec);
BfMat *bfMatBlockDenseToType(BfMat const *mat, BfType type);
BfSizeArray *bfMatBlockDenseGetNonzeroColumnRanges(BfMatBlockDense const *matBlockDense);
void bfMatBlockDensePrintBlocksDeep(BfMatBlockDense const *matBlockDense, FILE *fp, BfSize i0, BfSize j0, BfSize depth);

/** Interface: MatBlock */

BfSize bfMatBlockDenseNumBlocks(BfMatBlockDense const *matBlockDense);
BfSize bfMatBlockDenseGetNumRowBlocks(BfMatBlockDense const *matBlockDense);
BfSize bfMatBlockDenseGetNumColBlocks(BfMatBlockDense const *matBlockDense);
BfSize bfMatBlockDenseGetRowOffset(BfMatBlockDense const *matBlockDense, BfSize i);
BfSize bfMatBlockDenseGetColOffset(BfMatBlockDense const *matBlockDense, BfSize j);
BfMat const *bfMatBlockDenseGetBlockConst(BfMatBlockDense const *matBlockDense, BfSize i, BfSize j);

/** Implementation: MatBlockDense */

struct BfMatBlockDense {
  BfMatBlock super;
};

/* Upcasting: */
BfMat *bfMatBlockDenseToMat(BfMatBlockDense *matBlockDense);
BfMat const *bfMatBlockDenseConstToMatConst(BfMatBlockDense const *matBlockDense);
BfMatBlock *bfMatBlockDenseToMatBlock(BfMatBlockDense *matBlock);

/* Downcasting: */
BfMatBlockDense *bfMatToMatBlockDense(BfMat *mat);
BfMatBlockDense const *bfMatConstToMatBlockDenseConst(BfMat const *mat);

BfMatBlockDense *bfMatBlockDenseNew();
BfMatBlockDense *bfMatBlockDenseNewFromBlocks(BfSize numRowBlocks, BfSize numColBlocks, BfPtrArray *blocks);
BfMatBlockDense *bfMatBlockDenseNewRowFromBlocks(BfPtrArray *blocks);
BfMatBlockDense *bfMatBlockDenseNewColFromBlocks(BfPtrArray *blocks);
void bfMatBlockDenseInit(BfMatBlockDense *mat, BfSize numBlockRows,
                         BfSize numBlockCols);
void bfMatBlockDenseDeinit(BfMatBlockDense *mat);
void bfMatBlockDenseDealloc(BfMatBlockDense **mat);
void bfMatBlockDenseDeinitAndDealloc(BfMatBlockDense **mat);
BfMat *bfMatBlockDenseGetBlock(BfMatBlockDense *mat, BfSize i, BfSize j);
BfMat const *bfMatBlockDenseGetBlockConst(BfMatBlockDense const *mat, BfSize i, BfSize j);
void bfMatBlockDenseSetBlock(BfMatBlockDense *mat, BfSize i, BfSize j, BfMat *block);
void bfMatBlockDenseAppendColumn(BfMatBlockDense *mat, BfSize numBlockCols);
void bfMatBlockDenseAppendRow(BfMatBlockDense *mat, BfSize numBlockRows);
