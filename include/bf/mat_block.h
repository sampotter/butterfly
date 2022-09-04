#pragma once

#include "mat.h"

typedef struct BfMatBlockVtable {
  BfSize (*numBlocks)(BfMatBlock const *);
} BfMatBlockVtable;

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
   * `colOffset[numBlockCols]`, contains a sentinel value indicating
   * the total number of columns in the block matrix. */
  BfSize *colOffset;
};

BfMatBlock *bfMatToMatBlock(BfMat *mat);

void bfMatBlockInit(BfMatBlock *mat,
                    BfMatVtable *matVtbl, BfMatBlockVtable *matBlockVtbl,
                    BfSize numBlocks, BfSize numBlockRows, BfSize numBlockCols);
void bfMatBlockDeinit(BfMatBlock *mat);
void bfMatBlockDealloc(BfMatBlock **mat);
void bfMatBlockDeinitAndDealloc(BfMatBlock **mat);
BfSize bfMatBlockGetNumBlockRows(BfMatBlock const *mat, BfSize i);
BfSize bfMatBlockGetNumBlockCols(BfMatBlock const *mat, BfSize j);
BfSize bfMatBlockNumBlocks(BfMatBlock const *mat);
BfSize bfMatBlockGetNumRowBlocks(BfMatBlock const *mat);
BfSize bfMatBlockGetNumColBlocks(BfMatBlock const *mat);
BfSize bfMatBlockFindRowBlock(BfMatBlock const *mat, BfSize i);
BfSize bfMatBlockFindColBlock(BfMatBlock const *mat, BfSize i);

/* BfMat interface: */
BfMatType bfMatBlockGetType(BfMatBlock *mat);
BfSize bfMatBlockNumBytes(BfMatBlock *mat);
void bfMatBlockSave(BfMatBlock const *mat, char const *path);
BfMat *bfMatBlockMul(BfMatBlock const *op1, BfMat const *op2);
BfMat *bfMatBlockLstSq(BfMatBlock const *lhs, BfMat const *rhs);
