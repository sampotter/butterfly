#pragma once

#include "mat.h"

enum BfBlockMatProps {
  BF_BLOCK_MAT_PROP_NONE = 0,
  BF_BLOCK_MAT_PROP_DIAGONAL = (1 << 0),
  BF_BLOCK_MAT_PROP_VIEW = (1 << 1)
};

typedef struct BfBlockMat {
  enum BfDtypes dtype;
  enum BfBlockMatProps props;
  BfSize num_row_blocks, num_col_blocks;
  BfSize *row_sizes, *col_sizes;
  void *data;
} BfBlockMat;

enum BfError
bfInitEmptyBlockMat(BfBlockMat *block_mat,
                    enum BfDtypes dtype, enum BfBlockMatProps props,
                    BfSize num_row_blocks, BfSize *row_sizes,
                    BfSize num_col_blocks, BfSize *col_sizes);

BfMat *bfGetBlock(BfBlockMat *block, BfSize i, BfSize j);
