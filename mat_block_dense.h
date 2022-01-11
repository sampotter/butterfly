#pramga once

#include "mat_block.h"

typedef struct BfMatBlockDense {
  BfMatBlock super;
  BfMat **block;
} BfMatBlockDense;

BfMatBlockDense *bfMatBlockDenseNew();
void bfMatBlockInit(
void bfMatBlockDenseDeinit(BfMatBlockDense *mat);
void bfMatBlockDenseDelete(BfMatBlockDense **mat);
void bfMatBlockDenseDeinitAndDelete(BfMatBlockDense **mat);
BfMat *bfMatBlockDenseGetMatPtr(BfMatBlockDense *mat);
void bfMatBlockDenseSetBlock(BfMatBlockDense *mat, BfSize i, BfSize j, BfMat *block);
