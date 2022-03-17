#pragma once

#include "mat_block.h"

struct BfMatBlockDense {
  BfMatBlock super;
};

BF_DECLARE_INTERFACE_MAT(MatBlockDense);

BfMatBlockDense *bfMatBlockDenseNew();
void bfMatBlockDenseInit(BfMatBlockDense *mat, BfSize numBlockRows,
                         BfSize numBlockCols);
BfMat *bfMatBlockDenseGetMatPtr(BfMatBlockDense *mat);
BfMat const *bfMatBlockDenseGetMatConstPtr(BfMatBlockDense const *mat);
BfMat *bfMatBlockDenseGetBlock(BfMatBlockDense *mat, BfSize i, BfSize j);
void bfMatBlockDenseSetBlock(BfMatBlockDense *mat, BfSize i, BfSize j,
                             BfMat *block);

/* BfMatBlock interface */
BfSize bfMatBlockDenseNumBlocks(BfMatBlockDense const *mat);
