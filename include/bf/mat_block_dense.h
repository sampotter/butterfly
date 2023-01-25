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
void bfMatBlockDenseScaleCols(BfMat *mat, BfVec const *vec);
void bfMatBlockDenseAddInplace(BfMat *mat, BfMat const *otherMat);
BfMat *bfMatBlockDenseMul(BfMat const *mat, BfMat const *otherMat);
BfMat *bfMatBlockDenseToType(BfMat const *mat, BfType type);

/** Interface: MatBlock */

BfSize bfMatBlockDenseNumBlocks(BfMatBlock const *matBlock);
BfSize bfMatBlockDenseGetNumRowBlocks(BfMatBlock const *matBlock);
BfSize bfMatBlockDenseGetNumColBlocks(BfMatBlock const *matBlock);

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
