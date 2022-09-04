#pragma once

#include "mat.h"

#define BF_INTERFACE_MatBlock(Type, Subtype, _)                         \
  _(Type, Subtype, BfSize, NumBlocks, BfMatBlock const *)               \
  _(Type, Subtype, BfSize, GetNumRowBlocks, BfMatBlock const *)         \
  _(Type, Subtype, BfSize, GetNumColBlocks, BfMatBlock const *)         \
  _(Type, Subtype, BfSize, GetNumBlockRows, BfMatBlock const *, BfSize) \
  _(Type, Subtype, BfSize, GetNumBlockCols, BfMatBlock const *, BfSize)

#define INTERFACE BF_INTERFACE_MatBlock
BF_DEFINE_VTABLE_STRUCT(MatBlock);
BF_DECLARE_INTERFACE(MatBlock);
#undef INTERFACE

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

#define INTERFACE BF_INTERFACE_Mat
BF_DECLARE_INTERFACE(MatBlock);
#undef INTERFACE

#define INTERFACE BF_INTERFACE_MatBlock
BF_DECLARE_INTERFACE(MatBlock);
#undef INTERFACE

/* Upcasting: */
BfMat const *bfMatBlockConstToMatConst(BfMatBlock const *matBlock);

/* Downcasting: */
BfMatBlock *bfMatToMatBlock(BfMat *mat);
BfMatBlock const *bfMatConstToMatBlockConst(BfMat const *mat);

void bfMatBlockInit(BfMatBlock *mat,
                    BfMatVtable *matVtbl, BfMatBlockVtable *matBlockVtbl,
                    BfSize numBlocks, BfSize numBlockRows, BfSize numBlockCols);
void bfMatBlockDeinit(BfMatBlock *mat);
void bfMatBlockDealloc(BfMatBlock **mat);
void bfMatBlockDeinitAndDealloc(BfMatBlock **mat);
BfSize bfMatBlockFindRowBlock(BfMatBlock const *mat, BfSize i);
BfSize bfMatBlockFindColBlock(BfMatBlock const *mat, BfSize j);
