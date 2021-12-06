#include "block_mat.h"

#include <assert.h>
#include <stdlib.h>

/* this is messy... going to have to think about how to design a block
 * matrix type the right way */

enum BfError
bfInitEmptyBlockMat(BfBlockMat *block_mat,
                    enum BfDtypes dtype, enum BfBlockMatProps props,
                    BfSize num_row_blocks, BfSize *row_sizes,
                    BfSize num_col_blocks, BfSize *col_sizes)
{
  block_mat->dtype = dtype;
  block_mat->props = props;
  block_mat->num_row_blocks = num_row_blocks;
  block_mat->num_col_blocks = num_col_blocks;
  block_mat->row_sizes = row_sizes;
  block_mat->col_sizes = col_sizes;

  BfSize num_blocks;

  if (props & BF_MAT_PROP_DIAGONAL) {
    num_blocks = num_row_blocks<num_col_blocks?num_row_blocks:num_col_blocks;

    block_mat->data = malloc(num_blocks*sizeof(BfMat));

    BfMat *blocks = block_mat->data;
    for (BfSize i = 0; i < num_blocks; ++i) {
      bfInitEmptyMat(&blocks[i], dtype, BF_MAT_PROP_NONE,
                     (BfSize[]) {row_sizes[i], col_sizes[i]});
    }

  } else {
    return BF_ERROR_NOT_IMPLEMENTED;
  }

  return BF_ERROR_NO_ERROR;
}

BfMat *bfGetBlock(BfBlockMat *blockMat, BfSize i, BfSize j) {
  if (blockMat->props & BF_MAT_PROP_DIAGONAL) {
    if (i == j)
      return (BfMat *)blockMat->data + i;
  }

  assert(false);
}
