#pragma once

#include "mat_block.h"

/** Interface: Mat */

BfMat *bfMatBlockCooGetView(BfMatBlockCoo const *matBlockCoo);
BfMat *bfMatBlockCooCopy(BfMat const *mat);
BfMat *bfMatBlockCooSteal(BfMatBlockCoo *matBlockCoo);
BfVec *bfMatBlockCooGetRowCopy(BfMat const *mat, BfSize i);
void bfMatBlockCooDelete(BfMat **mat);
BfType bfMatBlockCooGetType(BfMat const *mat);
BfSize bfMatBlockCooNumBytes(BfMat const *mat);
void bfMatBlockCooDump(BfMatBlockCoo const *matBlockCoo, FILE *fp);
BfSize bfMatBlockCooGetNumRows(BfMat const *mat);
BfSize bfMatBlockCooGetNumCols(BfMat const *mat);
BfMat *bfMatBlockCooGetRowRangeCopy(BfMatBlockCoo const *matBlockCoo, BfSize i0, BfSize i1);
BfMat *bfMatBlockCooMul(BfMat const *op1, BfMat const *op2);
BfVec *bfMatBlockCooMulVec(BfMatBlockCoo const *matBlockCoo, BfVec const *vec);
BfVec *bfMatBlockCooRmulVec(BfMatBlockCoo const *matBlockCoo, BfVec const *vec);
void bfMatBlockCooNegate(BfMat *mat);
BfSizeArray *bfMatBlockCooGetNonzeroColumnRanges(BfMatBlockCoo const *matBlockCoo);
void bfMatBlockCooPrintBlocksDeep(BfMatBlockCoo const *matBlockCoo, FILE *fp, BfSize i0, BfSize j0, BfSize depth);

/** Interface: MatBlock */

BfSize bfMatBlockCooNumBlocks(BfMatBlockCoo const *matBlockCoo);
BfSize bfMatBlockCooGetNumRowBlocks(BfMatBlockCoo const *matBlockCoo);
BfSize bfMatBlockCooGetNumColBlocks(BfMatBlockCoo const *matBlockCoo);
BfSize bfMatBlockCooGetRowOffset(BfMatBlockCoo const *matBlockCoo, BfSize i);
BfSize bfMatBlockCooGetColOffset(BfMatBlockCoo const *matBlockCoo, BfSize j);
BfMat const *bfMatBlockCooGetBlockConst(BfMatBlockCoo const *matBlockCoo, BfSize i, BfSize j);

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

BfMatBlockCoo *bfMatBlockCooNew(void);
// BfMatBlockCoo *bfMatBlockCooNewFromArrays(BfSizeArray const *rowOffsets, BfSizeArray const *colOffsets, BfSizeArray const *rowInds, BfSizeArray const *colInds, BfPtrArray const *blocks);
BfMatBlockCoo *bfMatBlockCooNewColFromBlocks(BfPtrArray *blocks, BfPolicy policy);
BfMatBlockCoo *bfMatBlockCooNewRowFromBlocks(BfPtrArray *blocks, BfPolicy policy);
BfMatBlockCoo *bfMatBlockCooNewFromIndexedBlocks(BfSize numRows, BfSize numCols, BfPtrArray *indexedBlocks, BfPolicy policy);
void bfMatBlockCooInit(BfMatBlockCoo *mat, BfSize numBlockRows, BfSize numBlockCols, BfSize numBlocks);
void bfMatBlockCooDeinit(BfMatBlockCoo *mat);
void bfMatBlockCooDealloc(BfMatBlockCoo **mat);
void bfMatBlockCooDeinitAndDealloc(BfMatBlockCoo **mat);
