#pragma once

#include "mat.h"

/** Interface: MatBlock */

BfSize bfMatBlockNumBlocks(BfMatBlock const *);
BfSize bfMatBlockGetNumRowBlocks(BfMatBlock const *);
BfSize bfMatBlockGetNumColBlocks(BfMatBlock const *);
BfSize bfMatBlockGetNumBlockRows(BfMatBlock const *, BfSize);
BfSize bfMatBlockGetNumBlockCols(BfMatBlock const *, BfSize);
BfSize bfMatBlockGetRowOffset(BfMatBlock const *, BfSize);
BfSize bfMatBlockGetColOffset(BfMatBlock const *, BfSize);
BfMat *bfMatBlockGetBlock(BfMatBlock *, BfSize, BfSize);
BfMat const *bfMatBlockGetBlockConst(BfMatBlock const *, BfSize, BfSize);
BfMat *bfMatBlockGetBlockCopy(BfMatBlock const *, BfSize, BfSize);
void bfMatBlockGetRowSpan(BfMatBlock const *matBlock, BfSize i, BfSize *i0, BfSize *i1);
void bfMatBlockGetColSpan(BfMatBlock const *matBlock, BfSize j, BfSize *j0, BfSize *j1);

typedef struct BfMatBlockVtable {
  __typeof__(&bfMatBlockNumBlocks) NumBlocks;
  __typeof__(&bfMatBlockGetNumRowBlocks) GetNumRowBlocks;
  __typeof__(&bfMatBlockGetNumColBlocks) GetNumColBlocks;
  __typeof__(&bfMatBlockGetNumBlockRows) GetNumBlockRows;
  __typeof__(&bfMatBlockGetNumBlockCols) GetNumBlockCols;
  __typeof__(&bfMatBlockGetRowOffset) GetRowOffset;
  __typeof__(&bfMatBlockGetColOffset) GetColOffset;
  __typeof__(&bfMatBlockGetBlock) GetBlock;
  __typeof__(&bfMatBlockGetBlockConst) GetBlockConst;
  __typeof__(&bfMatBlockGetBlockCopy) GetBlockCopy;
  __typeof__(&bfMatBlockGetRowSpan) GetRowSpan;
  __typeof__(&bfMatBlockGetColSpan) GetColSpan;
} BfMatBlockVtable;

/** Implementation: MatBlock */

/* An abstract block matrix type. Shouldn't be instantiated
 * directly. */
struct BfMatBlock {
  BfMat super;

  BfMatBlockVtable *vtbl;

  /* The blocks. Their shapes must be compatible, but the type of each
   * block can vary. */
  BfMat **block;

  /* TODO: add the following sizes here, and store the real shape of
   * the matrix in BfMat. */
  // BfSize numRowBlocks;
  // BfSize numColBlocks;

  /* Array of `numBlockRows + 1` entries containing the offset in rows
   * of each block row in the matrix (i.e., `block[i]` starts on row
   * `rowOffset[rowInd[i]]`). The final entry,
   * `rowOffset[numBlockRows]`, contains a sentinel value indicating
   * the total number of rows in the block matrix. */
  BfSize *rowOffset;

  /* Array of `numBlockCols + 1` entries containing the offset in
   * columns of each block column in the matrix (i.e., `block[i]`
   * starts on column `colOffset[colInd[i]]`). The final entry,
   * `coloffset[numBlockCols]`, contains a sentinel value indicating
   * the total number of columns in the block matrix. */
  BfSize *colOffset;
};

/* Upcasting: */
BfMat const *bfMatBlockConstToMatConst(BfMatBlock const *matBlock);

/* Downcasting: */
BfMatBlock *bfMatToMatBlock(BfMat *mat);
BfMatBlock const *bfMatConstToMatBlockConst(BfMat const *mat);

void bfMatBlockInvalidate(BfMatBlock *matBlock);
void bfMatBlockInit(BfMatBlock *mat, BfMatVtable *matVtbl, BfMatBlockVtable *MAT_BLOCK_VTABLE, BfSize numBlocks, BfSize numBlockRows, BfSize numBlockCols);
void bfMatBlockDeinit(BfMatBlock *mat);
void bfMatBlockDealloc(BfMatBlock **mat);
void bfMatBlockDeinitAndDealloc(BfMatBlock **mat);
void bfMatBlockDelete(BfMat **mat);
BfSize bfMatBlockFindRowBlock(BfMatBlock const *mat, BfSize i);
BfSize bfMatBlockGetMaxDepth(BfMatBlock const *);
