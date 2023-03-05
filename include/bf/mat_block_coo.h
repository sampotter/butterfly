#pragma once

#include "mat_block.h"

/** Interface: Mat */

BfMat *bfMatBlockCooCopy(BfMat const *mat);
BfVec *bfMatBlockCooGetRowCopy(BfMat const *mat, BfSize i);
void bfMatBlockCooDelete(BfMat **mat);
BfType bfMatBlockCooGetType(BfMat const *mat);
BfSize bfMatBlockCooNumBytes(BfMat const *mat);
BfSize bfMatBlockCooGetNumRows(BfMat const *mat);
BfSize bfMatBlockCooGetNumCols(BfMat const *mat);
BfMat *bfMatBlockCooGetRowRangeCopy(BfMatBlockCoo const *matBlockCoo, BfSize i0, BfSize i1);
BfMat *bfMatBlockCooMul(BfMat const *op1, BfMat const *op2);
void bfMatBlockCooNegate(BfMat *mat);
void bfMatBlockCooPrintBlocksDeep(BfMatBlockCoo const *matBlockCoo, FILE *fp, BfSize i0, BfSize j0, BfSize depth);

/** Interface: MatBlock */

BfSize bfMatBlockCooNumBlocks(BfMatBlock const *mat);

struct BfMatBlockCoo {
  BfMatBlock super;

  BfSize numBlocks;

  /* The row which the block resides in (i.e., `block[i]` is in the
   * `rowInd[i]`th block row). */
  BfSize *rowInd;

  /* The column which the block resides in (i.e., `block[i]` in in the
   * `colInd[i]`th block column). */
  BfSize *colInd;
};

/** Upcasting: MatBlockCoo -> Mat */

BfMat *bfMatBlockCooToMat(BfMatBlockCoo *matBlockCoo);
BfMat const *bfMatBlockCooConstToMatConst(BfMatBlockCoo const *matBlockCoo);

/** Upcasting: MatBlockCoo -> MatBlock */

BfMatBlock *bfMatBlockCooToMatBlock(BfMatBlockCoo *matBlockCoo);
BfMatBlock const *bfMatBlockCooConstToMatBlockConst(BfMatBlockCoo const *matBlockCoo);

/** Downcasting: Mat -> MatBlockCoo */

BfMatBlockCoo *bfMatToMatBlockCoo(BfMat *mat);
BfMatBlockCoo const *bfMatConstToMatBlockCooConst(BfMat const *mat);

/** Downcasting: MatBlock -> MatBlockCoo */

BfMatBlockCoo const *bfMatBlockConstToMatBlockCooConst(BfMatBlock const *matBlock);

/** Implementation: MatBlockCoo */

BfMatBlockCoo *bfMatBlockCooNew();
BfMatBlockCoo *bfMatBlockCooNewFromArrays(BfSizeArray const *rowOffsets, BfSizeArray const *colOffsets, BfSizeArray const *rowInds, BfSizeArray const *colInds, BfPtrArray const *blocks);
BfMatBlockCoo *bfMatBlockCooNewColFromBlocks(BfPtrArray *blocks);
BfMatBlockCoo *bfMatBlockCooNewRowFromBlocks(BfPtrArray *blocks);
BfMatBlockCoo *bfMatBlockCooNewFromIndexedBlocks(BfSize numRows, BfSize numCols, BfPtrArray *indexedBlocks);
void bfMatBlockCooInit(BfMatBlockCoo *mat, BfSize numBlockRows, BfSize numBlockCols, BfSize numBlocks);
void bfMatBlockCooDeinit(BfMatBlockCoo *mat);
void bfMatBlockCooDealloc(BfMatBlockCoo **mat);
void bfMatBlockCooDeinitAndDealloc(BfMatBlockCoo **mat);
